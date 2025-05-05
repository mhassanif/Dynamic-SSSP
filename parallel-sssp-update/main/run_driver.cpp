// #include <mpi.h>
// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <vector>
// #include <tuple>
// #include <unordered_map>
// #include <set>
// #include <string>
// #include "../mpi/comm_utils.h"
// #include "../core/sssp_update.h"
// #include "../core/graph_structs.h"

// struct ChangeLog
// {
//     int global_id;
//     int old_distance;
//     int new_distance;
//     int old_parent;
//     int new_parent;
// };

// // Parse edge insertions and deletions from file
// void parse_updates(const std::string &insert_file, const std::string &delete_file,
//                    std::vector<std::tuple<int, int, int>> &insertions,
//                    std::vector<std::pair<int, int>> &deletions)
// {
//     std::ifstream ins(insert_file), dels(delete_file);
//     int u, v, w;

//     while (ins >> u >> v >> w)
//         insertions.emplace_back(u, v, w);

//     while (dels >> u >> v)
//         deletions.emplace_back(u, v);
// }

// int main(int argc, char **argv)
// {
//     MPI_Init(&argc, &argv);

//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     if (argc < 3)
//     {
//         if (rank == 0)
//             std::cerr << "Usage: mpirun -n <N> ./main <insertions.txt> <deletions.txt> [source=0]\n";
//         MPI_Finalize();
//         return 1;
//     }

//     std::string insert_file = argv[1];
//     std::string delete_file = argv[2];
//     int source = (argc > 3) ? std::stoi(argv[3]) : 0;

//     // ---------------------------
//     // Step 1: Load local partition
//     // ---------------------------
//     GraphPartition local_graph;
//     load_partition(rank, local_graph);

//     // Set source distance = 0 if local
//     auto src_it = local_graph.global_to_local.find(source);
//     if (src_it != local_graph.global_to_local.end())
//     {
//         int local_src = src_it->second;
//         local_graph.vertices[local_src].distance = 0;
//     }

//     // -----------------------------------------
//     // Step 2: Parse edge insertions/deletions
//     // -----------------------------------------
//     std::vector<std::tuple<int, int, int>> insertions;
//     std::vector<std::pair<int, int>> deletions;
//     parse_updates(insert_file, delete_file, insertions, deletions);

//     // ------------------------------------------------
//     // Step 3: Apply edge insertions to graph structure
//     // ------------------------------------------------
//     for (const auto &[u_gid, v_gid, weight] : insertions)
//     {
//         for (int gid : {u_gid, v_gid})
//         {
//             if (local_graph.global_to_local.count(gid) == 0)
//                 continue;
//             int local_idx = local_graph.global_to_local[gid];
//             Vertex &v = local_graph.vertices[local_idx];
//             int other_gid = (gid == u_gid) ? v_gid : u_gid;

//             // Avoid inserting duplicate edge (optional)
//             bool already_exists = false;
//             for (const Edge &e : v.edges)
//             {
//                 if (e.dest == other_gid)
//                 {
//                     already_exists = true;
//                     break;
//                 }
//             }
//             if (!already_exists)
//             {
//                 v.edges.push_back({other_gid, weight});
//             }
//         }
//     }

//     // ----------------------------------------
//     // Step 4: Save old distances/parents
//     // ----------------------------------------
//     std::unordered_map<int, std::pair<int, int>> old_state;
//     for (const auto &v : local_graph.vertices)
//         old_state[v.id] = {v.distance, v.parent};

//     // ----------------------------------------
//     // Step 5: Start timer and run update phase
//     // ----------------------------------------
//     double start_time = MPI_Wtime();

//     identifyAffectedVertices(local_graph, deletions, insertions);
//     propagateInfinity(local_graph);

//     bool converged = false;
//     while (!converged)
//     {
//         updateSSSP_OpenMP(local_graph, source);
//         exchange_boundary_data(local_graph);
//         converged = check_global_convergence(local_graph);

//         for (auto &v : local_graph.vertices)
//             v.updated = false;
//     }

//     double end_time = MPI_Wtime();

//     // ----------------------------------------
//     // Step 6: Collect change logs
//     // ----------------------------------------
//     std::vector<ChangeLog> local_changes;
//     for (const auto &v : local_graph.vertices)
//     {
//         auto it = old_state.find(v.id);
//         if (it != old_state.end())
//         {
//             int old_dist = it->second.first;
//             int old_par = it->second.second;
//             if (v.distance != old_dist || v.parent != old_par)
//             {
//                 local_changes.push_back({v.id, old_dist, v.distance, old_par, v.parent});
//             }
//         }
//     }

//     // ----------------------------------------
//     // Step 7: Gather logs at rank 0
//     // ----------------------------------------
//     std::vector<int> recv_counts(size);
//     int local_count = local_changes.size();
//     MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

//     std::vector<ChangeLog> all_changes;
//     std::vector<int> displs(size);
//     int total_changes = 0;

//     if (rank == 0)
//     {
//         displs[0] = 0;
//         for (int i = 1; i < size; ++i)
//             displs[i] = displs[i - 1] + recv_counts[i - 1];
//         total_changes = displs[size - 1] + recv_counts[size - 1];
//         all_changes.resize(total_changes);
//     }

//     MPI_Datatype MPI_ChangeLog;
//     MPI_Type_contiguous(5, MPI_INT, &MPI_ChangeLog);
//     MPI_Type_commit(&MPI_ChangeLog);

//     MPI_Gatherv(local_changes.data(), local_count, MPI_ChangeLog,
//                 all_changes.data(), recv_counts.data(), displs.data(),
//                 MPI_ChangeLog, 0, MPI_COMM_WORLD);

//     // ----------------------------------------
//     // Step 8: Print results
//     // ----------------------------------------
//     if (rank == 0)
//     {
//         std::cout << "-------------------------------------------\n";
//         std::cout << "SSSP Update Time: " << (end_time - start_time) << " seconds\n";
//         std::cout << "Affected Vertices: " << total_changes << "\n";
//         for (const auto &log : all_changes)
//         {
//             std::cout << "Vertex " << log.global_id
//                       << ": distance " << log.old_distance << " → " << log.new_distance
//                       << ", parent " << log.old_parent << " → " << log.new_parent << "\n";
//         }
//         std::cout << "-------------------------------------------\n";
//     }

//     MPI_Type_free(&MPI_ChangeLog);
//     MPI_Finalize();
//     return 0;
// }

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <chrono>
#include "../core/sssp_update.h"
#include "../core/graph_structs.h"
#include "../core/openmp_utils.h"

struct ChangeLog {
    int global_id;
    int old_distance;
    int new_distance;
    int old_parent;
    int new_parent;
};

void parse_updates(const std::string &insert_file, const std::string &delete_file,
                   std::vector<std::tuple<int, int, int>> &insertions,
                   std::vector<std::pair<int, int>> &deletions) {
    std::ifstream ins(insert_file), dels(delete_file);
    int u, v, w;

    while (ins >> u >> v >> w)
        insertions.emplace_back(u, v, w);

    while (dels >> u >> v)
        deletions.emplace_back(u, v);
}

void load_full_graph(const std::string &filename, GraphPartition &g) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Failed to open graph file: " << filename << std::endl;
        exit(1);
    }

    int num_vertices;
    infile >> num_vertices;
    g.initialize(num_vertices);
    g.local_to_global.resize(num_vertices);

    for (int i = 0; i < num_vertices; ++i) {
        int id, parent, dist;
        infile >> id >> parent >> dist;

        g.vertices[i].id = id;
        g.vertices[i].parent = parent;
        g.vertices[i].distance = dist;
        g.local_to_global[i] = id;
        g.global_to_local[id] = i;

        int neighbor, weight;
        while (infile.peek() != '\n' && infile >> neighbor >> weight) {
            g.vertices[i].edges.push_back({neighbor, weight});
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

int main(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "Usage: ./sssp_serial <graph_file.txt> <insertions.txt> <deletions.txt> [source=0]\n";
        return 1;
    }

    std::string graph_file = argv[1];
    std::string insert_file = argv[2];
    std::string delete_file = argv[3];
    int source = (argc > 4) ? std::stoi(argv[4]) : 0;

    GraphPartition g;
    load_full_graph(graph_file, g);

        // Dump graph structure for debugging
std::cout << "------ Initial Graph Structure ------\n";
for (const auto &v : g.vertices) {
    std::cout << "Vertex " << v.id << " (dist=" << v.distance << ", parent=" << v.parent << ") -> ";
    for (const auto &e : v.edges) {
        std::cout << "(" << e.dest << ", w=" << e.weight << ") ";
    }
    std::cout << "\n";
}
std::cout << "-------------------------------------\n";

    // Set source distance
    if (g.global_to_local.count(source))
        g.vertices[g.global_to_local[source]].distance = 0;

    // Load updates
    std::vector<std::tuple<int, int, int>> insertions;
    std::vector<std::pair<int, int>> deletions;
    parse_updates(insert_file, delete_file, insertions, deletions);

    // Save old state
    std::unordered_map<int, std::pair<int, int>> old_state;
    for (const auto &v : g.vertices)
        old_state[v.id] = {v.distance, v.parent};

    // Apply insertions
    for (const auto &[u, v, w] : insertions) {
        for (int a : {u, v}) {
            if (g.global_to_local.count(a)) {
                int idx = g.global_to_local[a];
                int b = (a == u) ? v : u;
                bool exists = false;
                for (auto &e : g.vertices[idx].edges)
                    if (e.dest == b) exists = true;
                if (!exists)
                    g.vertices[idx].edges.push_back({b, w});
            }
        }
    }

    // Apply deletions
    for (const auto &[u, v] : deletions) {
        for (int a : {u, v}) {
            if (g.global_to_local.count(a)) {
                int idx = g.global_to_local[a];
                int b = (a == u) ? v : u;
                auto &edges = g.vertices[idx].edges;
                edges.erase(std::remove_if(edges.begin(), edges.end(),
                    [&](Edge &e) { return e.dest == b; }), edges.end());
            }
        }
    }

    // Dump graph structure for debugging
std::cout << "------ post insertion/deletions Graph Structure ------\n";
for (const auto &v : g.vertices) {
    std::cout << "Vertex " << v.id << " (dist=" << v.distance << ", parent=" << v.parent << ") -> ";
    for (const auto &e : v.edges) {
        std::cout << "(" << e.dest << ", w=" << e.weight << ") ";
    }
    std::cout << "\n";
}
std::cout << "-------------------------------------\n";

    // Run update
    auto start = std::chrono::high_resolution_clock::now();

    identifyAffectedVertices(g, deletions, insertions);

    std::cout << "------ Affected Vertices After Identification ------\n";
for (const auto &v : g.vertices) {
    if (v.affected) {
        std::cout << "Vertex " << v.id << " marked as affected (distance: " 
                  << v.distance << ", parent: " << v.parent << ")\n";
    }
}
std::cout << "-----------------------------------------------------\n";



    propagateInfinity(g);

            // Dump graph structure for debugging
std::cout << "------ Initial Graph Structure ------\n";
for (const auto &v : g.vertices) {
    std::cout << "Vertex " << v.id << " (dist=" << v.distance << ", parent=" << v.parent << ") -> ";
    for (const auto &e : v.edges) {
        std::cout << "(" << e.dest << ", w=" << e.weight << ") ";
    }
    std::cout << "\n";
}
std::cout << "-------------------------------------\n";

    bool converged = false;
    while (!converged) {
        updateSSSP_OpenMP(g, source);
        converged = true;
        for (const auto &v : g.vertices) {
            if (v.updated) {
                converged = false;
                break;
            }
        }
        for (auto &v : g.vertices) v.updated = false;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();

    // Show result
    std::vector<ChangeLog> changes;
    for (const auto &v : g.vertices) {
        auto it = old_state.find(v.id);
        if (it != old_state.end()) {
            if (v.distance != it->second.first || v.parent != it->second.second) {
                changes.push_back({v.id, it->second.first, v.distance, it->second.second, v.parent});
            }
        }
    }

    std::cout << "-------------------------------------------\n";
    std::cout << "SSSP Update Time: " << duration << " seconds\n";
    std::cout << "Affected Vertices: " << changes.size() << "\n";
    for (const auto &log : changes) {
        std::cout << "Vertex " << log.global_id
                  << ": distance " << log.old_distance << " → " << log.new_distance
                  << ", parent " << log.old_parent << " → " << log.new_parent << "\n";
    }
    std::cout << "-------------------------------------------\n";

    return 0;
}
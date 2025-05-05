#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <set>
#include <string>
#include "../mpi/comm_utils.h"
#include "../core/sssp_update.h"
#include "../core/graph_structs.h"

struct ChangeLog
{
    int global_id;
    int old_distance;
    int new_distance;
    int old_parent;
    int new_parent;
};

// Parse edge insertions and deletions from file
void parse_updates(const std::string &insert_file, const std::string &delete_file,
                   std::vector<std::tuple<int, int, int>> &insertions,
                   std::vector<std::pair<int, int>> &deletions)
{
    std::ifstream ins(insert_file), dels(delete_file);
    int u, v, w;

    while (ins >> u >> v >> w)
        insertions.emplace_back(u, v, w);

    while (dels >> u >> v)
        deletions.emplace_back(u, v);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 3)
    {
        if (rank == 0)
            std::cerr << "Usage: mpirun -n <N> ./main <insertions.txt> <deletions.txt> [source=0]\n";
        MPI_Finalize();
        return 1;
    }

    std::string insert_file = argv[1];
    std::string delete_file = argv[2];
    int source = (argc > 3) ? std::stoi(argv[3]) : 0;

    // ---------------------------
    // Step 1: Load local partition
    // ---------------------------
    GraphPartition local_graph;
    load_partition(rank, local_graph);

    // Set source distance = 0 if local
    auto src_it = local_graph.global_to_local.find(source);
    if (src_it != local_graph.global_to_local.end())
    {
        int local_src = src_it->second;
        local_graph.vertices[local_src].distance = 0;
    }

    // -----------------------------------------
    // Step 2: Parse edge insertions/deletions
    // -----------------------------------------
    std::vector<std::tuple<int, int, int>> insertions;
    std::vector<std::pair<int, int>> deletions;
    parse_updates(insert_file, delete_file, insertions, deletions);

    // ------------------------------------------------
    // Step 3: Apply edge insertions to graph structure
    // ------------------------------------------------
    for (const auto &[u_gid, v_gid, weight] : insertions)
    {
        for (int gid : {u_gid, v_gid})
        {
            if (local_graph.global_to_local.count(gid) == 0)
                continue;
            int local_idx = local_graph.global_to_local[gid];
            Vertex &v = local_graph.vertices[local_idx];
            int other_gid = (gid == u_gid) ? v_gid : u_gid;

            // Avoid inserting duplicate edge (optional)
            bool already_exists = false;
            for (const Edge &e : v.edges)
            {
                if (e.dest == other_gid)
                {
                    already_exists = true;
                    break;
                }
            }
            if (!already_exists)
            {
                v.edges.push_back({other_gid, weight});
            }
        }
    }

    // ----------------------------------------
    // Step 4: Save old distances/parents
    // ----------------------------------------
    std::unordered_map<int, std::pair<int, int>> old_state;
    for (const auto &v : local_graph.vertices)
        old_state[v.id] = {v.distance, v.parent};

    // ----------------------------------------
    // Step 5: Start timer and run update phase
    // ----------------------------------------
    double start_time = MPI_Wtime();

    identifyAffectedVertices(local_graph, deletions, insertions);
    propagateInfinity(local_graph);

    bool converged = false;
    while (!converged)
    {
        updateSSSP_OpenMP(local_graph, source);
        exchange_boundary_data(local_graph);
        converged = check_global_convergence(local_graph);

        for (auto &v : local_graph.vertices)
            v.updated = false;
    }

    double end_time = MPI_Wtime();

    // ----------------------------------------
    // Step 6: Collect change logs
    // ----------------------------------------
    std::vector<ChangeLog> local_changes;
    for (const auto &v : local_graph.vertices)
    {
        auto it = old_state.find(v.id);
        if (it != old_state.end())
        {
            int old_dist = it->second.first;
            int old_par = it->second.second;
            if (v.distance != old_dist || v.parent != old_par)
            {
                local_changes.push_back({v.id, old_dist, v.distance, old_par, v.parent});
            }
        }
    }

    // ----------------------------------------
    // Step 7: Gather logs at rank 0
    // ----------------------------------------
    std::vector<int> recv_counts(size);
    int local_count = local_changes.size();
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<ChangeLog> all_changes;
    std::vector<int> displs(size);
    int total_changes = 0;

    if (rank == 0)
    {
        displs[0] = 0;
        for (int i = 1; i < size; ++i)
            displs[i] = displs[i - 1] + recv_counts[i - 1];
        total_changes = displs[size - 1] + recv_counts[size - 1];
        all_changes.resize(total_changes);
    }

    MPI_Datatype MPI_ChangeLog;
    MPI_Type_contiguous(5, MPI_INT, &MPI_ChangeLog);
    MPI_Type_commit(&MPI_ChangeLog);

    MPI_Gatherv(local_changes.data(), local_count, MPI_ChangeLog,
                all_changes.data(), recv_counts.data(), displs.data(),
                MPI_ChangeLog, 0, MPI_COMM_WORLD);

    // ----------------------------------------
    // Step 8: Print results
    // ----------------------------------------
    if (rank == 0)
    {
        std::cout << "-------------------------------------------\n";
        std::cout << "SSSP Update Time: " << (end_time - start_time) << " seconds\n";
        std::cout << "Affected Vertices: " << total_changes << "\n";
        for (const auto &log : all_changes)
        {
            std::cout << "Vertex " << log.global_id
                      << ": distance " << log.old_distance << " → " << log.new_distance
                      << ", parent " << log.old_parent << " → " << log.new_parent << "\n";
        }
        std::cout << "-------------------------------------------\n";
    }

    MPI_Type_free(&MPI_ChangeLog);
    MPI_Finalize();
    return 0;
}

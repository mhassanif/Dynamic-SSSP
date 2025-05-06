#include "comm_utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <mpi.h>

std::unordered_map<int, int> load_vertex_owner_map(const std::string &filepath)
{
    std::unordered_map<int, int> vertex_owner;
    std::ifstream infile(filepath);
    if (!infile)
    {
        std::cerr << "Failed to load vertex_owner file: " << filepath << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int owner;
    int global_id = 0;

    while (infile >> owner)
    {
        vertex_owner[global_id++] = owner;
    }

    return vertex_owner;
}
void load_partition(int rank, GraphPartition &g)
{
    // Load vertex owner map
    g.vertex_owner = load_vertex_owner_map("data/metis_output/division.txt");

    // Load partition data
    std::string filename = "data/metis_output/subgraph_" + std::to_string(rank) + ".txt";
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Rank " << rank << " failed to load file: " << filename << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int num_vertices;
    infile >> num_vertices;
    g.initialize(num_vertices);
    g.local_to_global.resize(num_vertices);

    // consume leftover newline
    std::string line;
    std::getline(infile, line);

    // Stage: read each line: global_id, global_parent, dist, [edges...] terminated by -1
    for (int i = 0; i < num_vertices; ++i) {
        if (!std::getline(infile, line)) {
            std::cerr << "Unexpected EOF reading vertex " << i << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::istringstream iss(line);
        int global_id, global_parent, dist;
        iss >> global_id >> global_parent >> dist;

        // assign id and distance
        g.vertices[i].id       = global_id;
        g.vertices[i].distance = dist;

        // convert global_parent → local index (or -1 if parent lives off‐rank)
        auto pit = g.global_to_local.find(global_parent);
        g.vertices[i].parent = (pit == g.global_to_local.end()
                                ? -1
                                : pit->second);

        // record mappings
        g.local_to_global[i]        = global_id;
        g.global_to_local[global_id] = i;

        // read edges until sentinel -1
        int dest, weight;
        while (iss >> dest && dest != -1) {
            if (!(iss >> weight)) {
                std::cerr << "Incomplete edge entry for vertex " << global_id << "\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            g.vertices[i].edges.push_back({ dest, weight });
        }
    }

    // Mark boundary vertices
    for (int i = 0; i < g.num_vertices; ++i) {
        auto &v = g.vertices[i];
        v.is_boundary = false;
        for (auto &e : v.edges) {
            int owner = g.vertex_owner[e.dest];
            if (owner != rank) {
                v.is_boundary = true;
                break;
            }
        }
    }
}



// Exchange distances of boundary vertices with other ranks
void exchange_boundary_data(GraphPartition &g)
{
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Prepare data to send: rank → (global_id, distance)
    std::unordered_map<int, std::vector<std::pair<int, int>>> outbound_data;

    for (const Vertex &v : g.vertices)
    {
        if (!v.is_boundary)
            continue;

        for (const Edge &e : v.edges)
        {
            int target_rank = g.vertex_owner[e.dest];
            if (target_rank != world_rank)
            {
                outbound_data[target_rank].emplace_back(v.id, v.distance);
            }
        }
    }

    // Exchange sizes
    std::vector<int> send_counts(world_size, 0);
    std::vector<int> recv_counts(world_size, 0);

    for (const auto &[rank, data] : outbound_data)
    {
        send_counts[rank] = static_cast<int>(data.size());
    }

    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Prepare send and receive buffers
    std::vector<std::vector<int>> send_buffers(world_size);
    std::vector<std::vector<int>> recv_buffers(world_size);

    for (int i = 0; i < world_size; ++i)
    {
        send_buffers[i].reserve(send_counts[i] * 2);
        recv_buffers[i].resize(recv_counts[i] * 2); // each entry is (global_id, distance)
    }

    for (int i = 0; i < world_size; ++i)
    {
        for (const auto &[gid, dist] : outbound_data[i])
        {
            send_buffers[i].push_back(gid);
            send_buffers[i].push_back(dist);
        }
    }

    // Perform MPI communication
    std::vector<MPI_Request> requests;

    for (int i = 0; i < world_size; ++i)
    {
        if (!send_buffers[i].empty())
        {
            MPI_Request req;
            MPI_Isend(send_buffers[i].data(), send_buffers[i].size(), MPI_INT, i, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
    }

    for (int i = 0; i < world_size; ++i)
    {
        if (!recv_buffers[i].empty())
        {
            MPI_Request req;
            MPI_Irecv(recv_buffers[i].data(), recv_buffers[i].size(), MPI_INT, i, 0, MPI_COMM_WORLD, &req);
            requests.push_back(req);
        }
    }

    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

    // Process received data to update ghost vertices
    for (int i = 0; i < world_size; ++i)
    {
        for (size_t j = 0; j < recv_buffers[i].size(); j += 2)
        {
            int gid = recv_buffers[i][j];
            int dist = recv_buffers[i][j + 1];

            auto it = g.global_to_local.find(gid);
            if (it != g.global_to_local.end())
            {
                int local_idx = it->second;
                if (local_idx >= 0 && local_idx < static_cast<int>(g.vertices.size()))
                {
                    Vertex &v = g.vertices[local_idx];
                    if (dist < v.distance)
                    {
                        v.distance = dist;
                        v.updated = true;
                    }
                }
            }
        }
    }
}

// Returns true if no updates occurred anywhere (i.e., converged)
bool check_global_convergence(const GraphPartition &g)
{
    int local_change = 0;

    for (const Vertex &v : g.vertices)
    {
        if (v.updated)
        {
            local_change = 1;
            break;
        }
    }

    int global_change = 0;
    MPI_Allreduce(&local_change, &global_change, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    return (global_change == 0); // True if no rank reported an update
}

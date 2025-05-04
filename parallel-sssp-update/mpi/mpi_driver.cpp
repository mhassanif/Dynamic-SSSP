#include <mpi.h>
#include <iostream>
#include "comm_utils.h"
#include "../core/sssp_update.h"
#include "../core/graph_structs.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    GraphPartition local_graph;
    load_partition(rank, local_graph);

    int source = 0;
    if (argc > 1)
    {
        source = std::stoi(argv[1]);
    }

    auto it = local_graph.global_to_local.find(source);
    if (it != local_graph.global_to_local.end())
    {
        int local_idx = it->second;
        local_graph.vertices[local_idx].distance = 0;
    }

    std::vector<int> changed_nodes; // for dynamic updates
    identifyAffectedVertices(local_graph, changed_nodes);
    propagateInfinity(local_graph);

    double start_time = MPI_Wtime();

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

    if (rank == 0)
    {
        std::cout << "SSSP converged in " << (end_time - start_time) << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}

#include <mpi.h>
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

    int source = 0; // Set source from CLI if needed

    // Initialize source vertex if owned by this rank
    int local_idx = local_graph.global_to_local[source];
    if (local_idx != -1)
    {
        local_graph.vertices[local_idx].distance = 0;
    }

    // Step 1: identify affected vertices (if this is a dynamic update)
    std::vector<int> changed_nodes; // empty or loaded based on use case
    identifyAffectedVertices(local_graph, changed_nodes);
    propagateInfinity(local_graph);

    // Main SSSP update loop
    bool converged = false;
    while (!converged)
    {
        updateSSSP_OpenMP(local_graph, source);

        exchange_boundary_data(local_graph); // ← must come before convergence check

        converged = check_global_convergence(local_graph);

        // Clear update flags for next iteration
        for (auto &v : local_graph.vertices)
            v.updated = false;
    }

    MPI_Finalize();
    return 0;
}

#include <mpi.h>
#include "comm_utils.h"
#include "../core/sssp_update.h"

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    GraphPartition local_graph;
    load_partition(rank, local_graph);

    int source = 0;                 // or from CLI
    std::vector<int> changed_nodes; // simulate or load

    identifyAffectedVertices(local_graph, changed_nodes);
    propagateInfinity(local_graph);

    bool converged = false;
    while (!converged)
    {
        updateSSSP_OpenMP(local_graph, source);
        converged = check_global_convergence(local_graph);
        exchange_boundary_data(local_graph);
    }

    MPI_Finalize();
    return 0;
}

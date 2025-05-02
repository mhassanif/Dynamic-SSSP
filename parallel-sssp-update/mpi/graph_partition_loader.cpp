#include "comm_utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>

void load_partition(int rank, GraphPartition &g)
{
    std::string filename = "partitions/partition_" + std::to_string(rank) + ".txt";
    std::ifstream infile(filename);
    if (!infile)
    {
        std::cerr << "Rank " << rank << " failed to load file: " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int num_vertices;
    infile >> num_vertices;
    g.initialize(num_vertices);

    for (int i = 0; i < num_vertices; ++i)
    {
        int u, dist, num_edges;
        infile >> u >> dist >> num_edges;
        g.vertices[i].id = u;
        g.vertices[i].distance = dist;
        for (int j = 0; j < num_edges; ++j)
        {
            int v, w;
            infile >> v >> w;
            g.vertices[i].edges.push_back({v, w});
        }
    }
}

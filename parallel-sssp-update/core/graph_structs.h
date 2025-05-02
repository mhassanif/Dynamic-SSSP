#ifndef GRAPH_STRUCTS_H
#define GRAPH_STRUCTS_H

#include <vector>
#include <limits>

struct Edge
{
    int dest;
    int weight;
};

struct Vertex
{
    int id;
    int distance;
    std::vector<Edge> edges;
    bool affected;
};

struct GraphPartition
{
    int num_vertices;
    std::vector<Vertex> vertices;
    std::vector<int> local_to_global;
    std::vector<int> global_to_local;

    void initialize(int n);
};

#endif

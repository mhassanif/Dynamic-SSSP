#ifndef GRAPH_STRUCTS_H
#define GRAPH_STRUCTS_H

#include <vector>
#include <unordered_map>
#include <limits>

struct Edge
{
    int dest;   // Global destination vertex ID
    int weight; // Edge weight
};

struct Vertex
{
    int id;                   // Global vertex ID
    int distance;             // Current shortest path estimate
    std::vector<Edge> edges;  // Outgoing edges
    bool affected = false;    // For dynamic SSSP updates
    bool updated = false;     // Flag: was distance updated this iteration?
    bool is_boundary = false; // True if it has neighbor in another partition
};

struct GraphPartition
{
    int num_vertices;
    std::vector<Vertex> vertices;
    std::vector<int> local_to_global;
    std::unordered_map<int, int> global_to_local;
    std::unordered_map<int, int> vertex_owner; // NEW: Global ID → owning rank

    void initialize(int n)
    {
        num_vertices = n;
        vertices.resize(n);
        local_to_global.resize(n);
        global_to_local.clear();
        vertex_owner.clear();
    }
};

#endif

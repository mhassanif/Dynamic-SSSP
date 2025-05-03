#ifndef CORE_GRAPH_STRUCTS_H
#define CORE_GRAPH_STRUCTS_H

#include <vector>
#include <unordered_map>
#include <limits>

static constexpr int INF_DIST = std::numeric_limits<int>::max();

// Edge in a partition, dest is the global vertex ID
struct Edge {
    int dest;   // Global destination vertex ID
    int weight; // Edge weight
};

// Vertex local structure
struct Vertex {
    int id;                   // Global vertex ID
    int parent;               // Local index of parent in SSSP tree (-1 if none)
    int distance;             // Current shortest path estimate
    std::vector<Edge> edges;  // Outgoing edges

    bool affected = false;    // Needs relaxation in update phase
    bool updated = false;     // Marked by deletion propagation
    bool is_boundary = false; // Boundary vertex flag
};

// Graph partition stored on each MPI rank
struct GraphPartition {
    int num_vertices;
    std::vector<Vertex> vertices;
    std::vector<int> local_to_global;            // local index -> global ID
    std::unordered_map<int,int> global_to_local; // global ID -> local index
    std::unordered_map<int,int> vertex_owner;    // global ID -> owning MPI rank

    // Initialize empty partition of size n
    void initialize(int n) {
        num_vertices = n;
        vertices.resize(n);
        local_to_global.resize(n);
        global_to_local.clear();
        vertex_owner.clear();

        for(int i = 0; i < n; ++i) {
            vertices[i].id = -1;
            vertices[i].parent = -1;
            vertices[i].distance = INF_DIST;
            vertices[i].edges.clear();
            vertices[i].affected = false;
            vertices[i].updated = false;
            vertices[i].is_boundary = false;
        }
    }
};

#endif // CORE_GRAPH_STRUCTS_H

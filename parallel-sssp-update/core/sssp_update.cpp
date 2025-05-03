#include "sssp_update.h"
#include "openmp_utils.h"
#include <algorithm>

void identifyAffectedVertices(
    GraphPartition& g,
    const std::vector<std::pair<int,int>>& deletions,
    const std::vector<std::tuple<int,int,int>>& insertions) {
    // Reset flags
    for(auto& v : g.vertices) {
        v.affected = false;
        v.updated = false;
    }

    // Process deletions (invalidate deeper endpoint)
    #pragma omp parallel for
    for(size_t i = 0; i < deletions.size(); ++i) {
        int u_gid = deletions[i].first;
        int v_gid = deletions[i].second;
        int lu = g.global_to_local[u_gid];
        int lv = g.global_to_local[v_gid];
        Vertex& Vu = g.vertices[lu];
        Vertex& Vv = g.vertices[lv];
        // If edge is part of tree
        if(Vv.parent == lu || Vu.parent == lv) {
            int y = (Vu.distance > Vv.distance ? lu : lv);
            g.vertices[y].distance = INF_DIST;
            g.vertices[y].parent = -1;
            g.vertices[y].updated = true;
            g.vertices[y].affected = true;
        }
    }

    // Process insertions (improve distances)
    #pragma omp parallel for
    for(size_t i = 0; i < insertions.size(); ++i) {
        int u_gid, v_gid, w;
        std::tie(u_gid, v_gid, w) = insertions[i];
        int lu = g.global_to_local[u_gid];
        int lv = g.global_to_local[v_gid];
        Vertex& Vu = g.vertices[lu];
        Vertex& Vv = g.vertices[lv];

        int x = (Vu.distance <= Vv.distance ? lu : lv);
        int y = (Vu.distance >  Vv.distance ? lu : lv);
        if(g.vertices[y].distance > g.vertices[x].distance + w) {
            g.vertices[y].distance = g.vertices[x].distance + w;
            g.vertices[y].parent = x;
            g.vertices[y].affected = true;
        }
    }
}

void propagateInfinity(GraphPartition& g) {
    bool again = true;
    while(again) {
        again = false;
        #pragma omp parallel for
        for(int i = 0; i < g.num_vertices; ++i) {
            Vertex& v = g.vertices[i];
            if(v.updated) {
                v.updated = false;
                for(const Edge& e : v.edges) {
                    int lc = g.global_to_local[e.dest];
                    Vertex& c = g.vertices[lc];
                    if(c.distance != INF_DIST) {
                        c.distance = INF_DIST;
                        c.parent = -1;
                        c.updated = true;
                        c.affected = true;
                        again = true;
                    }
                }
            }
        }
    }
}

void updateSSSP_OpenMP(GraphPartition& g, int source_gid) {
    // Initialize source in local partition
    int src = g.global_to_local[source_gid];
    g.vertices[src].distance = 0;
    g.vertices[src].parent = src;
    g.vertices[src].affected = true;

    bool again = true;
    while(again) {
        again = false;
        #pragma omp parallel for
        for(int i = 0; i < g.num_vertices; ++i) {
            Vertex& v = g.vertices[i];
            if(v.affected) {
                v.affected = false;
                for(const Edge& e : v.edges) {
                    int nl = g.global_to_local[e.dest];
                    Vertex& vn = g.vertices[nl];
                    int alt = v.distance + e.weight;
                    if(alt < vn.distance) {
                        vn.distance = alt;
                        vn.parent = i;
                        vn.affected = true;
                        again = true;
                    }
                }
            }
        }
    }
}

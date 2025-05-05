#include "sssp_update.h"
#include "openmp_utils.h"
#include <algorithm>

void identifyAffectedVertices(
    GraphPartition& g,
    const std::vector<std::pair<int, int>>& deletions,
    const std::vector<std::tuple<int, int, int>>& insertions)
{
    // Reset flags
    for (auto& v : g.vertices) {
        v.affected = false;
        v.updated = false;
    }

    // === Deletions (Undirected) ===
    #pragma omp parallel for
    for (size_t i = 0; i < deletions.size(); ++i) {
        int u_gid = deletions[i].first;
        int v_gid = deletions[i].second;

        if (g.global_to_local.count(u_gid) == 0 || g.global_to_local.count(v_gid) == 0)
            continue;  // Not in this partition

        int lu = g.global_to_local[u_gid];
        int lv = g.global_to_local[v_gid];
        Vertex& Vu = g.vertices[lu];
        Vertex& Vv = g.vertices[lv];

        // If edge was part of the tree, invalidate the deeper one
        if (Vu.parent == lv || Vv.parent == lu) {
            int y = (Vu.distance > Vv.distance ? lu : lv);
            g.vertices[y].distance = INF_DIST;
            g.vertices[y].parent = -1;
            g.vertices[y].updated = true;
            g.vertices[y].affected = true;
        }

        // Optional: remove the edge from both endpoints' edge lists
        // (You can skip this if deletions are logical only)
        Vu.edges.erase(std::remove_if(Vu.edges.begin(), Vu.edges.end(),
                        [v_gid](const Edge& e) { return e.dest == v_gid; }),
                        Vu.edges.end());
        Vv.edges.erase(std::remove_if(Vv.edges.begin(), Vv.edges.end(),
                        [u_gid](const Edge& e) { return e.dest == u_gid; }),
                        Vv.edges.end());
    }

    // === Insertions (Undirected) ===
    #pragma omp parallel for
    for (size_t i = 0; i < insertions.size(); ++i) {
        int u_gid, v_gid, w;
        std::tie(u_gid, v_gid, w) = insertions[i];

        if (g.global_to_local.count(u_gid) == 0 && g.global_to_local.count(v_gid) == 0)
            continue; // Neither endpoint is in this partition

        bool u_local = g.global_to_local.count(u_gid);
        bool v_local = g.global_to_local.count(v_gid);

        if (u_local) {
            int lu = g.global_to_local[u_gid];
            Vertex& Vu = g.vertices[lu];

            // Check if v_gid is reachable via mapping (remote or local)
            if (g.vertex_owner.count(v_gid)) {
                int lv = v_local ? g.global_to_local[v_gid] : -1;
                int v_dist = v_local ? g.vertices[lv].distance : INF_DIST;

                if (v_local && Vu.distance > v_dist + w) {
                    Vu.distance = v_dist + w;
                    Vu.parent = lv;
                    Vu.affected = true;
                }

                // Add the edge if it doesn't already exist
                bool exists = std::any_of(Vu.edges.begin(), Vu.edges.end(),
                    [v_gid](const Edge& e) { return e.dest == v_gid; });
                if (!exists) {
                    Vu.edges.push_back({v_gid, w});
                }
            }
        }

        if (v_local) {
            int lv = g.global_to_local[v_gid];
            Vertex& Vv = g.vertices[lv];

            if (g.vertex_owner.count(u_gid)) {
                int lu = u_local ? g.global_to_local[u_gid] : -1;
                int u_dist = u_local ? g.vertices[lu].distance : INF_DIST;

                if (u_local && Vv.distance > u_dist + w) {
                    Vv.distance = u_dist + w;
                    Vv.parent = lu;
                    Vv.affected = true;
                }

                // Add the edge if it doesn't already exist
                bool exists = std::any_of(Vv.edges.begin(), Vv.edges.end(),
                    [u_gid](const Edge& e) { return e.dest == u_gid; });
                if (!exists) {
                    Vv.edges.push_back({u_gid, w});
                }
            }
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

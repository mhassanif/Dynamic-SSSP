#include "sssp_update.h"
#include "openmp_utils.h"
#include <algorithm>
#include<iostream>

void identifyAffectedVertices(
    GraphPartition& g,
    const std::vector<std::pair<int, int>>& deletions,
    const std::vector<std::tuple<int, int, int>>& insertions)
{
    // Step 2: Reset all flags
    for (auto& v : g.vertices) {
        v.affected = false;
        v.updated = false;
        v.affectedDel = false;
    }

    // Step 3–9: Process deletions
    #pragma omp parallel for
    for (size_t i = 0; i < deletions.size(); ++i) {
        int u_gid = deletions[i].first;
        int v_gid = deletions[i].second;

        if (g.global_to_local.count(u_gid) == 0 || g.global_to_local.count(v_gid) == 0)
            continue;

        int lu = g.global_to_local[u_gid];
        int lv = g.global_to_local[v_gid];
        Vertex& Vu = g.vertices[lu];
        Vertex& Vv = g.vertices[lv];

        // If edge (u,v) is in T, remove deeper one
        if (Vu.parent == lv || Vv.parent == lu) {
            int y = (Vu.distance > Vv.distance ? lu : lv);
            g.vertices[y].distance = INF_DIST;
            g.vertices[y].parent = -1;
            g.vertices[y].updated = true;
            g.vertices[y].affected = true;
            g.vertices[y].affectedDel = true;
        }

        Vu.edges.erase(std::remove_if(Vu.edges.begin(), Vu.edges.end(),
                        [v_gid](const Edge& e) { return e.dest == v_gid; }),
                        Vu.edges.end());
        Vv.edges.erase(std::remove_if(Vv.edges.begin(), Vv.edges.end(),
                        [u_gid](const Edge& e) { return e.dest == u_gid; }),
                        Vv.edges.end());
    }

    // Step 10–19: Process insertions
    #pragma omp parallel for
    for (size_t i = 0; i < insertions.size(); ++i) {
        int u_gid, v_gid, w;
        std::tie(u_gid, v_gid, w) = insertions[i];

        bool u_local = g.global_to_local.count(u_gid);
        bool v_local = g.global_to_local.count(v_gid);

        int lu = u_local ? g.global_to_local[u_gid] : -1;
        int lv = v_local ? g.global_to_local[v_gid] : -1;

        int u_dist = u_local ? g.vertices[lu].distance : INF_DIST;
        int v_dist = v_local ? g.vertices[lv].distance : INF_DIST;

        int x_gid, y_gid, x_dist;
        if (u_dist <= v_dist) {
            x_gid = u_gid; y_gid = v_gid; x_dist = u_dist;
        } else {
            x_gid = v_gid; y_gid = u_gid; x_dist = v_dist;
        }

        bool y_local = g.global_to_local.count(y_gid);
        if (y_local && x_dist + w < g.vertices[g.global_to_local[y_gid]].distance) {
            int ly = g.global_to_local[y_gid];
            g.vertices[ly].distance = x_dist + w;
            g.vertices[ly].parent = g.global_to_local.count(x_gid) ? g.global_to_local[x_gid] : -1;
            g.vertices[ly].affected = true;
        }

        if (u_local) {
            auto& Vu = g.vertices[lu];
            if (std::none_of(Vu.edges.begin(), Vu.edges.end(), [v_gid](const Edge& e) { return e.dest == v_gid; }))
                Vu.edges.push_back({v_gid, w});
        }
        if (v_local) {
            auto& Vv = g.vertices[lv];
            if (std::none_of(Vv.edges.begin(), Vv.edges.end(), [u_gid](const Edge& e) { return e.dest == u_gid; }))
                Vv.edges.push_back({u_gid, w});
        }
    }
}




// Propagate infinite distance along tree-children per Algorithm 3 (lines 2–8)
// 2) propagateInfinity: guard every global_to_local lookup on e.dest
void propagateInfinity(GraphPartition& g) {
    // build children by local‐parent → local‐child
    std::vector<std::vector<int>> children(g.num_vertices);
    for (int i = 0; i < g.num_vertices; ++i) {
        int lp = g.vertices[i].parent;
        if (lp >= 0 && lp < g.num_vertices)
            children[lp].push_back(i);
    }

    bool again = true;
    while (again) {
        again = false;
        #pragma omp parallel for schedule(dynamic) reduction(||:again)
        for (int v = 0; v < g.num_vertices; ++v) {
            if (!g.vertices[v].affectedDel) continue;
            g.vertices[v].affectedDel = false;
            for (int c : children[v]) {
                auto &Cv = g.vertices[c];
                if (Cv.distance != INF_DIST) {
                    Cv.distance     = INF_DIST;
                    Cv.parent       = -1;
                    Cv.affected     = true;
                    Cv.affectedDel = true;
                    again = true;
                }
            }
        }
    }
}


// 3) updateSSSP_OpenMP: guard every e.dest lookup
void updateSSSP_OpenMP(GraphPartition& g, int source_gid) {
    // Safe source init
    auto it_src = g.global_to_local.find(source_gid);
    if (it_src != g.global_to_local.end()) {
        int src = it_src->second;
        g.vertices[src].distance     = 0;
        g.vertices[src].parent       = -1;
        g.vertices[src].affected     = true;
        g.vertices[src].affectedDel = false;
    }

    bool again = true;
    while (again) {
        again = false;
        #pragma omp parallel for schedule(dynamic) reduction(||:again)
        for (int i = 0; i < g.num_vertices; ++i) {
            auto &v = g.vertices[i];
            if (!v.affected) continue;
            v.affected = false;

            for (const auto &e : v.edges) {
                auto it_n = g.global_to_local.find(e.dest);
                if (it_n == g.global_to_local.end()) continue;
                int nl = it_n->second;
                auto &vn = g.vertices[nl];

                int alt = (v.distance == INF_DIST) ? INF_DIST : v.distance + e.weight;
                if (alt < vn.distance) {
                    vn.distance = alt;
                    vn.parent   = i;
                    vn.affected = true;
                    again       = true;
                }

                if (vn.distance != INF_DIST) {
                    int alt2 = vn.distance + e.weight;
                    if (alt2 < v.distance) {
                        v.distance = alt2;
                        v.parent   = nl;
                        v.affected = true;
                        again      = true;
                    }
                }
            }
        }
    }
}

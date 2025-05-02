#ifndef CORE_SSSP_UPDATE_H
#define CORE_SSSP_UPDATE_H

#include "graph_structs.h"
#include <vector>

// Algorithm 2: Mark vertices affected by deletions and insertions
void identifyAffectedVertices(
    GraphPartition& g,
    const std::vector<std::pair<int,int>>& deletions,
    const std::vector<std::tuple<int,int,int>>& insertions);

// Algorithm 3a: Propagate INF distances down deleted subtrees
void propagateInfinity(GraphPartition& g);

// Algorithm 3b: Relax all affected vertices in parallel (OpenMP)
void updateSSSP_OpenMP(GraphPartition& g, int source_gid);

#endif // CORE_SSSP_UPDATE_H

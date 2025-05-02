#ifndef SSSP_UPDATE_H
#define SSSP_UPDATE_H

#include "graph_structs.h"
#include <vector>

void identifyAffectedVertices(GraphPartition &g, const std::vector<int> &changed_vertices);
void propagateInfinity(GraphPartition &g);
void updateSSSP_OpenMP(GraphPartition &g, int source);

#endif

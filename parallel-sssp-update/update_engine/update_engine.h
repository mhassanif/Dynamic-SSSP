#ifndef UPDATE_ENGINE_H
#define UPDATE_ENGINE_H

#include "../core/graph_structs.h"
#include <iostream>
#include <string>
#include <vector>

void apply_edge_updates(GraphPartition &g, const std::vector<std::tuple<int, int, int>> &insertions,
                        const std::vector<std::pair<int, int>> &deletions);

#endif

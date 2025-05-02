#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include "../core/graph_structs.h"

void exchange_boundary_data(GraphPartition &g);
bool check_global_convergence(const GraphPartition &g);
void load_partition(int rank, GraphPartition &g);

#endif

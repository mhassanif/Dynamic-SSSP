#ifndef COMM_UTILS_H
#define COMM_UTILS_H

#include "../core/graph_structs.h"

// Load METIS-partitioned graph for a rank
void load_partition(int rank, GraphPartition &g);

// Exchange halo/boundary node distances
void exchange_boundary_data(GraphPartition &g);

// Check if any rank has changed distances (returns true if converged globally)
bool check_global_convergence(const GraphPartition &g);

#endif

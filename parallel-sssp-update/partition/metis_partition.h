#ifndef METIS_PARTITION_H
#define METIS_PARTITION_H

#include "../core/graph_structs.h"
#include <iostream>
#include <string>

void partition_graph_with_metis(GraphPartition &full_graph, int num_parts, const std::string &output_dir);

#endif

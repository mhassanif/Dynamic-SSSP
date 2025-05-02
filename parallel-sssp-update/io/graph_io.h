#ifndef GRAPH_IO_H
#define GRAPH_IO_H

#include <iostream>
#include <fstream>
#include <string>
#include "../core/graph_structs.h"

bool load_graph_from_file(const std::string &path, GraphPartition &g);
void save_graph_to_file(const std::string &path, const GraphPartition &g);

#endif

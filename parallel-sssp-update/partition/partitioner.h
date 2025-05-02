#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <metis.h>
#include <string>
#include <vector>

class Partitioner {
public:
    static void partitionGraph(const std::string& inputFile, const std::string& outputFile, idx_t numPartitions);
private:
    static void readGraph(const std::string& filename, idx_t& nvtxs, idx_t*& xadj, idx_t*& adjncy);
    static void writePartitions(const std::string& filename, idx_t* part, idx_t nvtxs);
};

#endif

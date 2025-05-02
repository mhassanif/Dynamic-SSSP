#include "partitioner.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input graph> <output partition> <num partitions>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    idx_t numPartitions = std::stoi(argv[3]);

    try {
        Partitioner::partitionGraph(inputFile, outputFile, numPartitions);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

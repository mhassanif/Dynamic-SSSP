#include "partitioner.h"
#include <iostream>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <num partitions>" << std::endl;
        return 1;
    }

    // Hardcoded input and output file paths
    std::string inputFile = "../data/metis_input/metis_graph.txt";
    std::string outputFile = "../data/metis_output/division.txt";

    idx_t numPartitions = std::stoi(argv[1]);

    try {
        Partitioner::partitionGraph(inputFile, outputFile, numPartitions);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

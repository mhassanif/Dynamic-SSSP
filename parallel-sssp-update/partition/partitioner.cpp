#include "partitioner.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

void Partitioner::partitionGraph(const std::string& inputFile, const std::string& outputFile, idx_t numPartitions) {
    idx_t nvtxs;      // Number of vertices
    idx_t* xadj;      // Adjacency index array
    idx_t* adjncy;    // Adjacency list array
    idx_t* adjwgt;    // Edge weight array
    idx_t* part;      // Partition vector
    idx_t edgecut;    // Number of edge cuts

    // Read the graph
    readGraph(inputFile, nvtxs, xadj, adjncy, adjwgt);

    // Allocate memory for partition vector
    part = new idx_t[nvtxs];

    // Prepare optional parameters
    idx_t ncon = 1;                 // Number of balancing constraints
    idx_t options[METIS_NOPTIONS];  // Options array
    METIS_SetDefaultOptions(options);

    // Partition the graph
    int result = METIS_PartGraphKway(
        &nvtxs, &ncon, xadj, adjncy, NULL, NULL, adjwgt, &numPartitions,
        NULL, NULL, options, &edgecut, part
    );

    if (result == METIS_OK) {
        std::cout << "Partitioning completed successfully. Edge cut: " << edgecut << std::endl;
        writePartitions(outputFile, part, nvtxs);
    } else {
        std::cerr << "Error during graph partitioning!" << std::endl;
    }

    // Free memory
    delete[] xadj;
    delete[] adjncy;
    delete[] adjwgt;
    delete[] part;
}

void Partitioner::readGraph(const std::string& filename, idx_t& nvtxs, idx_t*& xadj, idx_t*& adjncy, idx_t*& adjwgt) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open graph file.");
    }

    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    idx_t numEdges;

    iss >> nvtxs >> numEdges;

    xadj = new idx_t[nvtxs + 1];
    adjncy = new idx_t[2 * numEdges];
    adjwgt = new idx_t[2 * numEdges];

    idx_t edgeIndex = 0;
    xadj[0] = 0;

    for (idx_t i = 0; i < nvtxs; ++i) {
        std::getline(file, line);
        std::istringstream iss(line);
        idx_t neighbor, weight;

        while (iss >> neighbor >> weight) {
            adjncy[edgeIndex] = neighbor - 1; // Convert to 0-based index
            adjwgt[edgeIndex] = weight;
            ++edgeIndex;
        }

        xadj[i + 1] = edgeIndex;
    }

    file.close();
}

void Partitioner::writePartitions(const std::string& filename, idx_t* part, idx_t nvtxs) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open output file.");
    }

    for (idx_t i = 0; i < nvtxs; ++i) {
        file << part[i] << std::endl;
    }

    file.close();
}

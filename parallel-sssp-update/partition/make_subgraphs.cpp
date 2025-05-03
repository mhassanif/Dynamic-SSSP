#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>

using namespace std;

// Function to parse the input graph
void parseGraph(const string& filename, int& numVertices, vector<vector<pair<int, int>>>& graph) {
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    inFile >> numVertices;
    graph.resize(numVertices);

    int vertexId, distance, numEdges, neighbor, weight;
    for (int i = 0; i < numVertices; ++i) {
        inFile >> vertexId >> distance >> numEdges;
        for (int j = 0; j < numEdges; ++j) {
            inFile >> neighbor >> weight;
            graph[vertexId].emplace_back(neighbor, weight);
        }
    }

    inFile.close();
}

// Function to parse the partition vector
void parsePartition(const string& filename, vector<int>& partition) {
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    int part;
    while (inFile >> part) {
        partition.push_back(part);
    }

    inFile.close();
}

// Function to write a subgraph to a file
void writeSubgraph(const string& filename, const set<int>& vertices, const vector<vector<pair<int, int>>>& graph) {
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "Error opening file for writing: " << filename << endl;
        exit(1);
    }

    outFile << vertices.size() << endl;

    for (int vertex : vertices) {
        outFile << vertex << " 0 " << graph[vertex].size();
        for (const auto& edge : graph[vertex]) {
            if (vertices.count(edge.first)) {
                outFile << " " << edge.first << " " << edge.second;
            }
        }
        outFile << endl;
    }

    outFile.close();
}

// Main function to divide the graph into subgraphs
void divideGraph(const string& graphFile, const string& partitionFile, int numPartitions) {
    int numVertices;
    vector<vector<pair<int, int>>> graph;
    vector<int> partition;

    parseGraph(graphFile, numVertices, graph);
    parsePartition(partitionFile, partition);

    if (partition.size() != numVertices) {
        cerr << "Partition vector size does not match number of vertices." << endl;
        exit(1);
    }

    vector<set<int>> subgraphVertices(numPartitions);

    // Group vertices into subgraphs based on the partition vector
    for (int i = 0; i < numVertices; ++i) {
        subgraphVertices[partition[i]].insert(i);
    }

    // Write each subgraph to a file
    for (int i = 0; i < numPartitions; ++i) {
        string filename = "subgraph_" + to_string(i) + ".txt";
        writeSubgraph(filename, subgraphVertices[i], graph);
        cout << "Subgraph " << i << " written to " << filename << endl;
    }
}

int main() {
    string graphFile = "../data/preprocessed/graph1.txt";          // Input graph file
    string partitionFile = "../data/metis_output/division.txt"; // Partition vector file
    int numPartitions = 4;                  // Number of partitions (adjust as needed)

    divideGraph(graphFile, partitionFile, numPartitions);

    return 0;
}

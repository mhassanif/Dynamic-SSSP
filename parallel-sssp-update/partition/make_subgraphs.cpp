#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <sstream>

using namespace std;

// Function to parse the input graph
void parseGraph(const string& filename, int& numVertices, vector<vector<string>>& graph) {
    ifstream inFile(filename);
    if (!inFile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    string line;
    getline(inFile, line); // First line contains numVertices and numEdges
    stringstream ss(line);
    ss >> numVertices;

    graph.resize(numVertices);

    int vertexId = 0;
    while (getline(inFile, line)) {
        graph[vertexId].push_back(line); // Store the entire line as a single string
        ++vertexId;
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
void writeSubgraph(const string& filename, const set<int>& vertices, const vector<vector<string>>& graph) {
    ofstream outFile(filename);
    if (!outFile) {
        cerr << "Error opening file for writing: " << filename << endl;
        exit(1);
    }

    outFile << vertices.size() << " " << "placeholder_edges_count" << endl; // Placeholder for edge count
    for (int vertex : vertices) {
        outFile << graph[vertex][0] << endl; // Write the full line corresponding to the vertex
    }

    outFile.close();
}

// Main function to divide the graph into subgraphs
void divideGraph(const string& graphFile, const string& partitionFile, int numPartitions) {
    int numVertices;
    vector<vector<string>> graph;
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

    // Write each subgraph to a file in the metis_output folder
    for (int i = 0; i < numPartitions; ++i) {
        string filename = "../data/metis_output/subgraph_" + to_string(i) + ".txt";
        writeSubgraph(filename, subgraphVertices[i], graph);
        cout << "Subgraph " << i << " written to " << filename << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <num partitions>" << endl;
        return 1;
    }

    int numPartitions = stoi(argv[1]);

    string graphFile = "../data/preprocessed/graph1.txt";          // Input graph file
    string partitionFile = "../data/metis_output/division.txt"; // Partition vector file

    divideGraph(graphFile, partitionFile, numPartitions);

    return 0;
}

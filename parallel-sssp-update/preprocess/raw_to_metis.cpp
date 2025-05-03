#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>

using namespace std;

int main() {
    string inputFile = "../data/raw/graph1.txt"; // Input raw graph file
    string outputFile = "../data/metis_input/metis_graph.txt"; // Output METIS graph file
    ifstream in(inputFile);
    ofstream out(outputFile);

    if (!in || !out) {
        cerr << "Error opening input or output file.\n";
        return 1;
    }

    int numVertices, numEdges;
    in >> numVertices >> numEdges;
    in.ignore(); // Move to next line

    // Adjacency list with weights
    vector<map<int, int>> adj(numVertices);

    // Read and populate the adjacency list
    for (int u = 0; u < numVertices; ++u) {
        string line;
        getline(in, line);
        istringstream iss(line);
        int v, w;
        while (iss >> v >> w) {
            // Only insert u->v if not already stored to avoid duplicate undirected edges
            if (adj[v].find(u) == adj[v].end()) {
                adj[u][v] = w;
            }
        }
    }

    // Count unique edges (each undirected edge once)
    int trueEdges = 0;
    for (int u = 0; u < numVertices; ++u) {
        trueEdges += adj[u].size();
    }

    // Write METIS header: <#vertices> <#edges> <format_flag>
    out << numVertices << " " << trueEdges << " 1\n";

    // Write adjacency list (1-based indices)
    for (int u = 0; u < numVertices; ++u) {
        for (const auto& [v, w] : adj[u]) {
            out << v + 1 << " " << w << " ";
        }
        out << "\n";
    }

    in.close();
    out.close();

    cout << "Conversion complete. Output written to: " << outputFile << "\n";
    return 0;
}
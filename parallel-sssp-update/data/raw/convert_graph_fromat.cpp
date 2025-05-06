#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <string>

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: ./convert_bitcoinotc_to_graph_format <input.csv> <output.txt>\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    // Step 1: Read the dataset and collect unique vertices and edges
    std::ifstream infile(input_file);
    if (!infile) {
        std::cerr << "Failed to open input file: " << input_file << "\n";
        return 1;
    }

    std::set<int> unique_vertices;
    std::vector<std::tuple<int, int, int>> edges; // (source, target, adjusted_weight)
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(iss, token, ','))
            tokens.push_back(token);

        if (tokens.size() != 4) {
            std::cerr << "Invalid line format: " << line << "\n";
            continue;
        }

        int source = std::stoi(tokens[0]);
        int target = std::stoi(tokens[1]);
        int rating = std::stoi(tokens[2]);
        // Adjust weight: RATING + 11 to ensure no zero weights (maps [-10, +10] to [1, 21])
        int weight = rating + 11;

        unique_vertices.insert(source);
        unique_vertices.insert(target);
        edges.emplace_back(source, target, weight);
    }
    infile.close();

    int num_vertices = unique_vertices.size();
    int num_edges = edges.size(); // Directed edges

    // Step 2: Map original vertex IDs to contiguous range [0, num_vertices-1]
    std::unordered_map<int, int> old_to_new_id;
    int new_id = 0;
    for (int v : unique_vertices) {
        old_to_new_id[v] = new_id++;
    }

    // Step 3: Build adjacency list for each vertex
    std::vector<std::vector<std::pair<int, int>>> adj_list(num_vertices);
    for (const auto& [source, target, weight] : edges) {
        int new_source = old_to_new_id[source];
        int new_target = old_to_new_id[target];
        adj_list[new_source].emplace_back(new_target, weight);
    }

    // Step 4: Write the graph to the output file
    std::ofstream outfile(output_file);
    if (!outfile) {
        std::cerr << "Failed to open output file: " << output_file << "\n";
        return 1;
    }

    // Write header: num_vertices num_edges
    outfile << num_vertices << " " << num_edges << "\n";

    // Write each vertex's adjacency list
    for (int v = 0; v < num_vertices; ++v) {
        const auto& neighbors = adj_list[v];
        if (neighbors.empty()) {
            outfile << "\n"; // Empty line for vertices with no outgoing edges
            continue;
        }
        for (size_t i = 0; i < neighbors.size(); ++i) {
            outfile << neighbors[i].first << " " << neighbors[i].second;
            if (i < neighbors.size() - 1)
                outfile << " ";
        }
        outfile << "\n";
    }

    outfile.close();
    std::cout << "Conversion complete. Output written to " << output_file << "\n";
    return 0;
}
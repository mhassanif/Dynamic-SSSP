#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <limits>
#include <string>
#include <chrono> // For timing

const int INF = std::numeric_limits<int>::max();

struct Edge
{
    int to;
    int weight;
};

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: ./dijkstra-on-raw <input_file> <output_file>\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    std::ifstream infile(input_file);
    if (!infile)
    {
        std::cerr << "Failed to open input file: " << input_file << "\n";
        return 1;
    }

    int num_vertices, num_edges;
    infile >> num_vertices >> num_edges;

    std::vector<std::vector<Edge>> adj(num_vertices);

    std::string line;
    std::getline(infile, line); // Consume leftover newline

    for (int u = 0; u < num_vertices; ++u)
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        int v, w;
        while (iss >> v >> w)
        {
            adj[u].push_back({v, w});
        }
    }
    infile.close();

    std::vector<int> dist(num_vertices, INF);
    std::vector<int> parent(num_vertices, -1);
    dist[0] = 0;

    using pii = std::pair<int, int>; // (distance, vertex)
    std::priority_queue<pii, std::vector<pii>, std::greater<>> pq;
    pq.push({0, 0});

    // Start timing
    auto start_time = std::chrono::high_resolution_clock::now();

    while (!pq.empty())
    {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u])
            continue;

        for (const auto &edge : adj[u])
        {
            int v = edge.to, w = edge.weight;
            if (dist[u] + w < dist[v])
            {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.push({dist[v], v});
            }
        }
    }

    // End timing
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Dijkstra execution time: " << elapsed.count() << " seconds\n";

    std::ofstream outfile(output_file);
    if (!outfile)
    {
        std::cerr << "Failed to open output file: " << output_file << "\n";
        return 1;
    }

    outfile << num_vertices << " " << num_edges << "\n";
    for (int u = 0; u < num_vertices; ++u)
    {
        outfile << u << " " << parent[u] << " ";
        outfile << (dist[u] == INF ? -1 : dist[u]);

        for (const auto &edge : adj[u])
        {
            outfile << " " << edge.to << " " << edge.weight;
        }
        outfile << "\n";
    }

    outfile.close();
    std::cout << "Dijkstra preprocessing complete. Output written to " << output_file << "\n";
    return 0;
}
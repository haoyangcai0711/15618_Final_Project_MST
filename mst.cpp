#include <vector>
#include <queue>
#include <iostream>
#include <climits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cstdlib>
#include <chrono>

struct CSRGraph {
    std::vector<int> values;
    std::vector<int> col_indices;
    std::vector<int> row_ptr;
};

struct Edge {
    int u;
    int v;
    int weight;
    bool operator<(const Edge& other) const {
        return weight < other.weight;
    }
};

std::vector<Edge> primMST(const CSRGraph& graph) {
    int numVertices = graph.row_ptr.size() - 1;
    std::vector<int> parent(numVertices, -1);
    std::vector<int> key(numVertices, INT_MAX);
    std::vector<bool> inMST(numVertices, false);
    std::vector<Edge> mstEdges;

    std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

    int startVertex = 0;
    key[startVertex] = 0;
    pq.push(std::make_pair(0, startVertex));

    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();

        if (inMST[u])
            continue;

        inMST[u] = true;

        if (parent[u] != -1) {
            Edge edge = {parent[u], u, key[u]};
            mstEdges.push_back(edge);
            // mstEdges.push_back(std::make_pair(parent[u], u));
        }

        for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
            int v = graph.col_indices[i];
            int weight = graph.values[i];

            if (!inMST[v] && weight < key[v]) {
                parent[v] = u;
                key[v] = weight;
                pq.push(std::make_pair(key[v], v));
            }
        }
    }

    return mstEdges;
}

struct UnionFind {
    std::vector<int> parent;
    std::vector<int> rank;

    UnionFind(int n) : parent(n), rank(n, 0) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
        }
    }

    int find(int u) {
        if (u != parent[u]) {
            parent[u] = find(parent[u]); // Path compression
        }
        return parent[u];
    }

    void unionSets(int u, int v) {
        int rootU = find(u);
        int rootV = find(v);
        if (rootU != rootV) {
            if (rank[rootU] > rank[rootV]) {
                parent[rootV] = rootU;
            } else if (rank[rootU] < rank[rootV]) {
                parent[rootU] = rootV;
            } else {
                parent[rootV] = rootU;
                rank[rootU]++;
            }
        }
    }
};

std::vector<Edge> extractEdgesFromCSR(const CSRGraph& graph) {
    std::vector<Edge> edges;
    int num_vertices = graph.row_ptr.size() - 1;
    for (int u = 0; u < num_vertices; ++u) {
        for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
            int v = graph.col_indices[i];
            int weight = graph.values[i];
            if (u < v) { // To avoid duplicate edges in undirected graph
                edges.push_back({u, v, weight});
            }
        }
    }
    return edges;
}

std::vector<Edge> kruskalMST(const CSRGraph& graph) {
    std::vector<Edge> edges = extractEdgesFromCSR(graph);
    std::sort(edges.begin(), edges.end());
    UnionFind uf(graph.row_ptr.size() - 1);
    std::vector<Edge> result;

    for (const Edge& e : edges) {
        if (uf.find(e.u) != uf.find(e.v)) {
            result.push_back(e);
            uf.unionSets(e.u, e.v);
        }
    }

    return result;
}

std::vector<int> read_array_from_file(std::ifstream& file, int size) {
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::vector<int> array;
    array.reserve(size);
    int num;
    while (iss >> num) {
        array.push_back(num);
    }
    return array;
}

int main(int argc, char* argv[]) {
    int scale = 3; // can be 1 - 5, the graph contains 10^scale vertices
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "--scale" || arg == "-s") && i + 1 < argc) {
            scale = std::atoi(argv[i + 1]);
            i++;
        }
    }
    std::string filename = "graph_" + std::to_string(scale) + ".txt";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::vector<int> graph_size = read_array_from_file(file, 2);
    int numEdges = graph_size[0];
    int numVertices = graph_size[1];

    std::vector<int> data = read_array_from_file(file, numEdges);
    std::vector<int> col_indices = read_array_from_file(file, numEdges);
    std::vector<int> row_ptr = read_array_from_file(file, numVertices + 1);

    file.close();

    CSRGraph graph;
    graph.values = data;
    graph.col_indices = col_indices;
    graph.row_ptr = row_ptr;
    
    auto startKruskal = std::chrono::high_resolution_clock::now();
    std::vector<Edge> mstEdgesKruskal = kruskalMST(graph);
    auto endKruskal = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedKruskal = endKruskal - startKruskal;
    std::cout << "Kruskal's algorithm took " << elapsedKruskal.count() << " milliseconds.\n";

    auto startPrim = std::chrono::high_resolution_clock::now();
    std::vector<Edge> mstEdgesPrim = primMST(graph);
    auto endPrim = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedPrim = endPrim - startPrim;
    std::cout << "Prim's algorithm took " << elapsedPrim.count() << " milliseconds.\n";
    
    int totalWeightKruskal = 0;
    for (const auto& edge : mstEdgesKruskal) {
        for (int i = graph.row_ptr[edge.u]; i < graph.row_ptr[edge.u + 1]; ++i) {
            if (graph.col_indices[i] == edge.v) {
                totalWeightKruskal += graph.values[i];
                break;
            }
        }
    }
    std::cout << "Total weight of MST (Kruskal): " << totalWeightKruskal << std::endl;

    int totalWeightPrim = 0;
    for (const auto& edge : mstEdgesPrim) {
        for (int i = graph.row_ptr[edge.u]; i < graph.row_ptr[edge.u + 1]; ++i) {
            if (graph.col_indices[i] == edge.v) {
                totalWeightPrim += graph.values[i];
                break;
            }
        }
    }
    std::cout << "Total weight of MST (Prim): " << totalWeightPrim << std::endl;

    return 0;
}


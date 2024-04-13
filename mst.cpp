#include <vector>
#include <queue>
#include <iostream>
#include <climits>
#include <fstream>
#include <sstream>
#include <algorithm>
//hi
struct CSRGraph {
    std::vector<int> values;
    std::vector<int> col_indices;
    std::vector<int> row_ptr;
};

std::vector<std::pair<int, int>> primMST(const CSRGraph& graph) {
    int numVertices = graph.row_ptr.size() - 1;
    std::vector<int> parent(numVertices, -1);
    std::vector<int> key(numVertices, INT_MAX);
    std::vector<bool> inMST(numVertices, false);
    std::vector<std::pair<int, int>> mstEdges;

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
            mstEdges.push_back(std::make_pair(parent[u], u));
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

// void primMST(const CSRGraph& graph, std::vector<std::pair<int, int>>& mstEdges) {
//     int num_vertices = graph.row_ptr.size() - 1;
//     std::vector<bool> inMST(num_vertices, false);
//     std::vector<int> key(num_vertices, INT_MAX);
//     std::vector<int> parent(num_vertices, -1);  // To store the MST
//     std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<std::pair<int, int>>> pq;

//     // Start with the first vertex in MST.
//     key[0] = 0;
//     pq.push({0, 0});  // weight, vertex

//     while (!pq.empty()) {
//         int u = pq.top().second;
//         pq.pop();

//         if (inMST[u]) continue;
//         inMST[u] = true;

//         // Iterate over all adjacent vertices of u
//         for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
//             int v = graph.col_indices[i];
//             int weight = graph.values[i];

//             // If v is not in MST and weight of (u,v) is smaller
//             if (!inMST[v] && key[v] > weight) {
//                 parent[v] = u;
//                 key[v] = weight;
//                 pq.push({key[v], v});
//             }
//         }
//     }

//     mstEdges.clear();
//     for (int v = 1; v < num_vertices; ++v) {
//         if (parent[v] != -1) {
//             mstEdges.push_back({parent[v], v});
//         }
//     }

//     std::cout << "Edges in the MST:" << std::endl;
//     for (const auto& edge : mstEdges) {
//         std::cout << edge.first << " - " << edge.second << std::endl;
//     }
// }

class DisjointSet {
private:
    std::vector<int> parent;
    std::vector<int> rank;

public:
    DisjointSet(int size) {
        parent.resize(size);
        rank.resize(size, 0);
        for (int i = 0; i < size; ++i) {
            parent[i] = i;
        }
    }

    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    void unionSet(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX != rootY) {
            if (rank[rootX] < rank[rootY]) {
                parent[rootX] = rootY;
            } else if (rank[rootX] > rank[rootY]) {
                parent[rootY] = rootX;
            } else {
                parent[rootY] = rootX;
                rank[rootX]++;
            }
        }
    }
};

std::vector<std::pair<int, int>> kruskalMST(const CSRGraph& graph) {
    int numVertices = graph.row_ptr.size() - 1;
    std::vector<std::pair<int, std::pair<int, int>>> edges;

    // Extract edges from the CSR graph
    for (int u = 0; u < numVertices; ++u) {
        for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
            int v = graph.col_indices[i];
            int weight = graph.values[i];
            edges.push_back(std::make_pair(weight, std::make_pair(u, v)));
        }
    }

    // Sort edges in ascending order of weights
    std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

    std::vector<std::pair<int, int>> mstEdges;
    DisjointSet ds(numVertices);

    for (const auto& edge : edges) {
        int u = edge.second.first;
        int v = edge.second.second;

        if (ds.find(u) != ds.find(v)) {
            ds.unionSet(u, v);
            mstEdges.push_back(std::make_pair(u, v));
        }
    }

    return mstEdges;
}

// struct Edge {
//     int u;
//     int v;
//     int weight;
//     bool operator<(const Edge& other) const {
//         return weight < other.weight;
//     }
// };

// struct UnionFind {
//     std::vector<int> parent;
//     std::vector<int> rank;

//     UnionFind(int n) : parent(n), rank(n, 0) {
//         for (int i = 0; i < n; ++i) {
//             parent[i] = i;
//         }
//     }

//     int find(int u) {
//         if (u != parent[u]) {
//             parent[u] = find(parent[u]); // Path compression
//         }
//         return parent[u];
//     }

//     void unionSets(int u, int v) {
//         int rootU = find(u);
//         int rootV = find(v);
//         if (rootU != rootV) {
//             if (rank[rootU] > rank[rootV]) {
//                 parent[rootV] = rootU;
//             } else if (rank[rootU] < rank[rootV]) {
//                 parent[rootU] = rootV;
//             } else {
//                 parent[rootV] = rootU;
//                 rank[rootU]++;
//             }
//         }
//     }
// };

// std::vector<Edge> extractEdgesFromCSR(const CSRGraph& graph) {
//     std::vector<Edge> edges;
//     int num_vertices = graph.row_ptr.size() - 1;
//     for (int u = 0; u < num_vertices; ++u) {
//         for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
//             int v = graph.col_indices[i];
//             int weight = graph.values[i];
//             if (u < v) { // To avoid duplicate edges in undirected graph
//                 edges.push_back({u, v, weight});
//             }
//         }
//     }
//     return edges;
// }

// std::vector<Edge> kruskalMST(const CSRGraph& graph) {
//     std::vector<Edge> edges = extractEdgesFromCSR(graph);
//     std::sort(edges.begin(), edges.end());
//     UnionFind uf(graph.row_ptr.size() - 1);
//     std::vector<Edge> result;

//     for (const Edge& e : edges) {
//         if (uf.find(e.u) != uf.find(e.v)) {
//             result.push_back(e);
//             uf.unionSets(e.u, e.v);
//         }
//     }

//     return result;
// }

std::vector<int> read_array_from_file(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    std::istringstream iss(line);
    std::vector<int> array;
    int num;
    while (iss >> num) {
        array.push_back(num);
    }
    return array;
}

int main() {
    std::ifstream file("csr_graph.txt");
    if (!file.is_open()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    std::vector<int> data = read_array_from_file(file);
    std::vector<int> col_indices = read_array_from_file(file);
    std::vector<int> row_ptr = read_array_from_file(file);

    file.close();

    CSRGraph graph;
    graph.values = data;
    graph.col_indices = col_indices;
    graph.row_ptr = row_ptr;

    std::vector<std::pair<int, int>> mstEdgesPrim = primMST(graph);
    std::vector<std::pair<int, int>> mstEdgesKruskal = kruskalMST(graph);

    std::cout << "Edges in the MST (Prim):" << std::endl;
    for (const auto& edge : mstEdgesPrim) {
        std::cout << edge.first << " - " << edge.second << std::endl;
    }
    std::cout << "Edges in the MST (Kruskal):" << std::endl;
    for (const auto& edge : mstEdgesKruskal) {
        std::cout << edge.first << " - " << edge.second << std::endl;
    }

    int totalWeightPrim = 0;
    for (const auto& edge : mstEdgesPrim) {
        totalWeightPrim += graph.values[graph.row_ptr[edge.first] + edge.second];
    }
    int totalWeightKruskal = 0;
    for (const auto& edge : mstEdgesKruskal) {
        totalWeightKruskal += graph.values[graph.row_ptr[edge.first] + edge.second];
    }
    std::cout << "Total weight of MST (Prim): " << totalWeightPrim << std::endl;
    std::cout << "Total weight of MST (Kruskal): " << totalWeightKruskal << std::endl;
}

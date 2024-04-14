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
#include <mutex>
#include <omp.h>
#include <atomic>
#include <cassert>

struct CSRGraph {
    std::vector<int> values; 
    std::vector<int> col_indices;
    std::vector<int> row_ptr; 
};

struct Edge {
    int u;
    int v;
    int weight;
    Edge(int u, int v, int weight) : u(u), v(v), weight(weight) {}
    // ~Edge() {
    //     std::cout << "Edge destructor called: " << u << " " << v << " " << weight << std::endl;
    // }
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

    bool unionSets(int u, int v) {
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
            return true;
        }
        return false;
    }
};

struct UnionFindThreadSafe {
    std::vector<int> parent;
    std::vector<int> rank;
    std::atomic_flag lock = ATOMIC_FLAG_INIT;

    bool tryLock() {
        return !lock.test_and_set(std::memory_order_acquire);
    }

    void unlock() {
        lock.clear(std::memory_order_release);
    }

    UnionFindThreadSafe(int n) : parent(n), rank(n) {
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    int find(int u) {
        // Path compression using iterative approach
        int root = u;
        while (root != parent[root]) {
            root = parent[root];
        }
        while (u != root) {
            int next = parent[u];
            parent[u] = root;
            u = next;
        }
        return root;
    }

    bool unionSets(int u, int v) {
        while (true) {
            u = find(u);
            v = find(v);

            if (u == v) {
                return false;
            }
            
            if (rank[u] > rank[v]) {
                if (!tryLock()) continue;
                parent[v] = u;
                unlock();
                return true;
            } else if (rank[u] < rank[v]) {
                if (!tryLock()) continue;
                parent[u] = v;
                unlock();
                return true;
            } else {
                if (!tryLock()) continue;
                parent[v] = u;
                if (rank[u] == rank[v]) {
                    rank[u]++;
                }
                unlock();
                return true;
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

std::vector<Edge> boruvkaMST(const CSRGraph& graph) {
    int numVertices = graph.row_ptr.size() - 1;
    UnionFind uf(numVertices);
    std::vector<Edge> mstEdges;
    int numComponents = numVertices;

    while (numComponents > 1) {
        std::vector<Edge> cheapest(numVertices, { -1, -1, std::numeric_limits<int>::max() });
        // Find the cheapest outgoing edge for each component
        for (int u = 0; u < numVertices; ++u) {
            int compU = uf.find(u);
            for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
                int v = graph.col_indices[i];
                int weight = graph.values[i];
                int compV = uf.find(v);
                if (compU != compV && weight < cheapest[compU].weight) {
                    cheapest[compU] = { u, v, weight };
                }
            }
        }
        // Add the cheapest outgoing edges to the MST
        for (int i = 0; i < numVertices; ++i) {
            Edge edge = cheapest[i];
            if (edge.u != -1 && uf.find(edge.u) != uf.find(edge.v)) {
                mstEdges.push_back(edge);
                uf.unionSets(edge.u, edge.v);
                numComponents--;
            }
        }
    }
    return mstEdges;
}

std::vector<Edge> boruvkaMSTParallel(const CSRGraph& graph, int numThreads = 8, int granularity = 8) {
    int numVertices = graph.row_ptr.size() - 1;
    UnionFindThreadSafe uf(numVertices);
    std::vector<Edge> mstEdges;
    int numComponents = numVertices;

    while (numComponents > 1) {
        std::vector<Edge> cheapest(numVertices, { -1, -1, std::numeric_limits<int>::max() });
        // Find the cheapest outgoing edge for each component
        int u = 0;
        #pragma omp parallel for default(none) private(u) shared(graph, uf, cheapest, numVertices) schedule(dynamic, 8) num_threads(numThreads)
        for (u = 0; u < numVertices; ++u) {
            int compU = uf.find(u);
            for (int i = graph.row_ptr[u]; i < graph.row_ptr[u + 1]; ++i) {
                int v = graph.col_indices[i];
                int weight = graph.values[i];
                int compV = uf.find(v);
                if (compU != compV && weight < cheapest[compU].weight) {
                    cheapest[compU] = { u, v, weight };
                }
            }
        }
        // Add the cheapest outgoing edges to the MST
        for (int i = 0; i < numVertices; ++i) {
            Edge edge = cheapest[i];
            if (edge.u != -1 && uf.find(edge.u) != uf.find(edge.v)) {
                mstEdges.push_back(edge);
                uf.unionSets(edge.u, edge.v);
                numComponents--;
            }
        }
    }
    return mstEdges;
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

// Usage: ./mst -s <graph size, exponent of 10> -t <number of threads>
int main(int argc, char* argv[]) {
    int scale = 3; // can be 1 - 5, the graph contains 10^scale vertices
    int numThreads = 8;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if ((arg == "--scale" || arg == "-s") && i + 1 < argc) {
            scale = std::atoi(argv[i + 1]);
            i++;
        } else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
            numThreads = std::atoi(argv[i + 1]);
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

    auto startBoruvka = std::chrono::high_resolution_clock::now();
    std::vector<Edge> mstEdgesBoruvka = boruvkaMST(graph);
    auto endBoruvka = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedBoruvka = endBoruvka - startBoruvka;
    std::cout << "Boruvka's algorithm took " << elapsedBoruvka.count() << " milliseconds.\n";

    auto startBoruvkaParallel = std::chrono::high_resolution_clock::now();
    std::vector<Edge> mstEdgesBoruvkaParallel = boruvkaMSTParallel(graph, numThreads);
    auto endBoruvkaParallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedBoruvkaParallel = endBoruvkaParallel - startBoruvkaParallel;
    std::cout << "Boruvka's algorithm (parallel) took " << elapsedBoruvkaParallel.count() << " milliseconds.\n" << std::flush;
    
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

    int totalWeightBoruvka = 0;
    for (const auto& edge : mstEdgesBoruvka) {
        for (int i = graph.row_ptr[edge.u]; i < graph.row_ptr[edge.u + 1]; ++i) {
            if (graph.col_indices[i] == edge.v) {
                totalWeightBoruvka += graph.values[i];
                break;
            }
        }
    }
    std::cout << "Total weight of MST (Boruvka): " << totalWeightBoruvka << std::endl;

    int totalWeightBoruvkaParallel = 0;
    for (const auto& edge : mstEdgesBoruvkaParallel) {
        for (int i = graph.row_ptr[edge.u]; i < graph.row_ptr[edge.u + 1]; ++i) {
            if (graph.col_indices[i] == edge.v) {
                totalWeightBoruvkaParallel += graph.values[i];
                break;
            }
        }
    }
    std::cout << "Total weight of MST (Boruvka parallel): " << totalWeightBoruvkaParallel << std::endl;

    return 0;
}


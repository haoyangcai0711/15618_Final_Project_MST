#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <queue>
#include <shared_mutex>
#include <sstream>
#include <string>
#include <vector>

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
  //     std::cout << "Edge destructor called: " << u << " " << v << " " <<
  //     weight << std::endl;
  // }
  bool operator<(const Edge &other) const { return weight < other.weight; }
};

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

std::vector<Edge> kruskalMST(std::vector<Edge> &edges, int numNodes) {
  std::sort(edges.begin(), edges.end());
  UnionFind uf(numNodes);
  std::vector<Edge> result;

  for (const Edge &e : edges) {
    if (uf.find(e.u) != uf.find(e.v)) {
      result.push_back(e);
      uf.unionSets(e.u, e.v);
    }
  }

  return result;
}

std::vector<Edge> extractEdgesFromCSR(const CSRGraph &graph) {
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

int totalWeightMST(const std::vector<std::vector<Edge>> &mst) {
  int total_weight = 0;
  for (const std::vector<Edge> &edges : mst) {
    for (const Edge &e : edges) {
      total_weight += e.weight;
    }
  }
  return total_weight / 2;
}

bool dfs(std::vector<std::vector<Edge>> &mst, int u, int v,
         std::vector<bool> &visited, std::deque<int> &path) {
  visited[u] = true;
  path.push_back(u);
  if (u == v) {
    return true;
  }
  for (const Edge &e : mst[u]) {
    if (!visited[e.v]) {
      if (dfs(mst, e.v, v, visited, path)) {
        return true;
      }
    }
  }
  path.pop_back();
  return false;
}

bool bfs(std::vector<std::vector<Edge>> &mst,
         std::vector<std::shared_mutex> &mst_locks, int u, int v,
         std::vector<bool> &visited, std::deque<Edge> &path) {
  std::queue<int> q;
  std::vector<Edge> parent(mst.size(), {INT_MAX, INT_MAX, INT_MAX});

  visited[u] = true;
  q.push(u);

  bool found = false;
  while (!q.empty() && !found) {
    int current = q.front();
    q.pop();
    std::shared_lock lock(mst_locks[current]);
    for (Edge &e : mst[current]) {
      if (!visited[e.v]) {
        visited[e.v] = true;
        parent[e.v] = e;
        q.push(e.v);
        if (e.v == v) {
          found = true;
          break;
        }
      }
    }
  }

  if (found) {
    // Reconstruct path from v to u using the parent pointers
    int at = v;
    while (at != u) {
      Edge edge = parent[at];
      path.push_front(edge);
      at = edge.u; // Move to the next node in the path
    }
  }

  return found;
}

void updateMST(std::vector<std::vector<Edge>> &mst,
               std::vector<std::shared_mutex> &mst_locks, Edge &new_edge) {
  // int original_weight = totalWeightMST(mst);
  while (true) {

    std::vector<bool> visited(mst.size(), false);

    std::deque<Edge> path_edges;
    if (!bfs(mst, mst_locks, new_edge.u, new_edge.v, visited, path_edges)) {
      printf("Thread %d: No path found between %d and %d (The original path "
             "had been "
             "modified by other threads)\n",
             omp_get_thread_num(), new_edge.u, new_edge.v);
      continue;
    }

    Edge max_edge = path_edges[0];
    for (auto e : path_edges) {
      if (e.weight > max_edge.weight) {
        max_edge = e;
      }
    }

    if (max_edge.weight > new_edge.weight) {
      std::unique_lock lock1(mst_locks[max_edge.u]);
      auto it = std::find_if(
          mst[max_edge.u].begin(), mst[max_edge.u].end(),
          [&max_edge](const Edge &e) { return e.v == max_edge.v; });
      if (it == mst[max_edge.u].end()) {
        printf(
            "Thread %d: Max Edge %d-%d %d has been removed by another thread\n",
            omp_get_thread_num(), max_edge.u, max_edge.v, max_edge.weight);
        continue;
      }

      std::unique_lock lock2(mst_locks[max_edge.v]);
      auto it2 = std::find_if(
          mst[max_edge.v].begin(), mst[max_edge.v].end(),
          [&max_edge](const Edge &e) { return e.v == max_edge.u; });
      if (it2 == mst[max_edge.v].end()) {
        printf(
            "Thread %d: Max Edge %d-%d %d has been removed by another thread\n",
            omp_get_thread_num(), max_edge.v, max_edge.u, max_edge.weight);
        continue;
      }

      mst[max_edge.u].erase(it);
      mst[max_edge.v].erase(it2);
      lock1.unlock();
      lock2.unlock();

      std::unique_lock lock3(mst_locks[new_edge.u]);
      std::unique_lock lock4(mst_locks[new_edge.v]);
      mst[new_edge.u].push_back(new_edge);
      mst[new_edge.v].push_back({new_edge.v, new_edge.u, new_edge.weight});
      lock3.unlock();
      lock4.unlock();
    }
    break;
  }
}

std::vector<int> read_array_from_file(std::ifstream &file, int size) {
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

int main(int argc, char *argv[]) {
  int scale = 3; // can be 1 - 5, the graph contains 10^scale vertices
  // int numThreads = 8;
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if ((arg == "--scale" || arg == "-s") && i + 1 < argc) {
      scale = std::atoi(argv[i + 1]);
      i++;
    }
    // else if ((arg == "--threads" || arg == "-t") && i + 1 < argc) {
    //     numThreads = std::atoi(argv[i + 1]);
    //     i++;
    // }
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

  std::vector<Edge> edges = extractEdgesFromCSR(graph);
  std::random_shuffle(edges.begin(), edges.end());
  // split the edges into 2 parts, with ratio 0.9 vs 0.1.
  int split = edges.size() * 0.95;
  std::vector<Edge> edges_95(edges.begin(), edges.begin() + split);
  std::vector<Edge> edges_05(edges.begin() + split, edges.end());

  auto startKruskal = std::chrono::high_resolution_clock::now();
  std::vector<Edge> mstEdgesKruskal = kruskalMST(edges_95, numVertices);
  auto endKruskal = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsedKruskal =
      endKruskal - startKruskal;
  std::cout << "Kruskal's algorithm took " << elapsedKruskal.count()
            << " milliseconds.\n";

  std::vector<std::vector<Edge>> mst(numVertices);
  for (Edge e : mstEdgesKruskal) {
    mst[e.u].push_back(e);
    mst[e.v].push_back({e.v, e.u, e.weight});
  }

  std::vector<std::shared_mutex> mst_locks(mst.size());
  auto startUpdate = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
  for (Edge &e : edges_05) {
    updateMST(mst, mst_locks, e);
  }
  auto endUpdate = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> elapsedUpdate =
      endUpdate - startUpdate;
  std::cout << "Updating MST took " << elapsedUpdate.count()
            << " milliseconds.\n";

  // verify the correctness of the updated MST
  int total_weight = totalWeightMST(mst);
  std::cout << "Total weight of the updated MST: " << total_weight << std::endl;

  // recompute the MST using Kruskal's algorithm
  std::vector<Edge> mstEdgesKruskal2 = kruskalMST(edges, numVertices);
  int total_weight_kruskal = 0;
  for (const Edge &e : mstEdgesKruskal2) {
    total_weight_kruskal += e.weight;
  }
  std::cout << "Total weight of the MST computed by Kruskal's algorithm: "
            << total_weight_kruskal << std::endl;

  return 0;
}

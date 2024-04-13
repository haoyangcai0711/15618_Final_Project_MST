# MST Algorithms

This repository contains tools and implementations for exploring and testing different Minimum Spanning Tree (MST) algorithms. It includes a Python script to generate random connected Watts Strogatz graphs and a C++ program to compute MSTs using Kruskal's and Prim's algorithms.

## Prerequisites

Ensure you have Python 3 and a C++ compiler installed on your system. This project also requires `make` for building the C++ code.

## Generating Graphs

The script `generate_graph.py` is used to generate random connected Watts Strogatz graphs. Here's how you can use it:

### Usage

```bash
python3 generate_graph.py -n <number_of_nodes> -k <nearest_neighbors> -p <rewiring_prob> -f <output_file>
```

- `-n`: Number of nodes in the graph.
- `-k`: Each node is joined with its k nearest neighbors in a ring topology.
- `-p`: The probability of rewiring each edge.
- `-f`: Output file name, which follows the naming convention graph_i.txt, where i indicates that the graph contains 10^i nodes.

### Example

```bash
python3 generate_graph.py -n 100000 -k 100 -p 0.3 -f graph_5.txt
```

## Building and Running the MST Algorithms

To compile the MST implementation, use the included Makefile:
```bash
make
```

After building the project, you can run the MST algorithms using:
```
./mst -s <size>
```
- `-s`: Specifies the size of the graph as an exponent of 10. For example, -s 5 means the graph size is 10^5 and the program will read from graph_5.txt.
This command will execute both Kruskal's and Prim's algorithms on the graph described in graph_\<size\>.txt, outputting the runtimes and total weights of the MSTs generated (the total weights should match).

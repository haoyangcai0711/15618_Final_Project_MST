import networkx as nx
import numpy as np
import scipy.sparse as sp

def generate_random_connected_graph(n, k, p, visualize=False, file_path="graph.txt", image_path="graph.png"):
    # Generate a connected Watts-Strogatz small-world graph
    G = nx.connected_watts_strogatz_graph(n, k, p)
    for (u, v, w) in G.edges(data=True):
        w['weight'] = np.random.randint(1, 50)
    
    # Save the graph to a file
    sparse_matrix = nx.to_scipy_sparse_matrix(G, format='csr')
    row_ptr = sparse_matrix.indptr
    col_ind = sparse_matrix.indices
    data = sparse_matrix.data
    with open(file_path, 'w') as f:
        f.write(f"{G.number_of_edges()} {G.number_of_nodes()}\n")
        f.write(" ".join(map(str, row_ptr)) + "\n")
        f.write(" ".join(map(str, col_ind)) + "\n")
        f.write(" ".join(map(str, data)) + "\n")
    
    # Compute the MST
    mst = nx.minimum_spanning_tree(G, algorithm='prim')
    total_weight = sum(data['weight'] for (u, v, data) in mst.edges(data=True))
    print(f"Graph saved to {file_path}, total weight of MST: {total_weight}")
    
    if visualize:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(16, 8))
        # Plot the original graph
        plt.subplot(121)
        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, node_color='orange', edge_color='gray', node_size=700)
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
        plt.title('Random Connected Graph')
        # Plot the MST
        plt.subplot(122)
        nx.draw(mst, pos, with_labels=True, node_color='lightblue', edge_color='blue', node_size=700)
        mst_edge_labels = nx.get_edge_attributes(mst, 'weight')
        nx.draw_networkx_edge_labels(mst, pos, edge_labels=mst_edge_labels)
        plt.title('Minimum Spanning Tree')
        plt.savefig(image_path)
        plt.savefig(image_path)

# Read arguments from command line
import argparse
parser = argparse.ArgumentParser(description='Generate a random connected graph')
parser.add_argument('-n', type=int, default=10, help='Number of nodes')
parser.add_argument('-k', type=int, default=4, help='Each node is connected to k nearest neighbors in ring topology')
parser.add_argument('-p', type=float, default=0.3, help='Probability of rewiring each edge')
parser.add_argument('-v', action='store_true', help='Visualize the graph')
parser.add_argument('-f', type=str, default='graph.txt', help='Output file path')
parser.add_argument('-i', type=str, default='graph.png', help='Output image path')

args = parser.parse_args()
generate_random_connected_graph(args.n, args.k, args.p, args.v, args.f, args.i)
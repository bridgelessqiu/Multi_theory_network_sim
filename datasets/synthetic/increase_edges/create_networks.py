# Erdös-Rényi networks

import networkx as nx
import numpy as np

# Number of vertices
n = 5000

# list of p
p_list = [0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25]

for p in p_list:
    # Graph
    G = nx.gnp_random_graph(n, p)

    # The path
    path = 'er_' + str(p) + '.edges'

    # Save the edgelist file
    nx.write_edgelist(G, path, delimiter = ' ', data = False)
# Erdös-Rényi networks

import networkx as nx
import numpy as np

# Edge creation probability for ER graphs
p = 0.1

low = 1000
high = 50000
step = 3000

for n in range(low, high, step):
    # Graph
    G = nx.gnp_random_graph(n,p)

    # The path
    path = 'er_' + str(n) + '.edges'

    # Save the edgelist file
    nx.write_edgelist(G, path, delimiter = ' ', data = False)
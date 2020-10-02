import networkx as nx
import numpy as np

"""
A: n x n adjacency matrix
T: n x c threshold matrix: between [0, d(v)+2]
B: n x c configuration matrix: random from 0 or 1
W: c x c weight matrix
d: c x 1 density vector: between [1, c]
"""

noise_levels = ['0', '0.01', '0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25']

network_name = 'google'

for noise in noise_levels:
    # file path
    f_path = network_name + '_' + noise + '_g2.edges'

    # the number of contagions
    c = 3

    # read in the graph
    G = nx.read_edgelist(f_path, nodetype = int, data=False)

    # the number of vertices
    n = G.number_of_nodes()

    ####################
    # Adjacency matrix #
    ####################
    A = nx.to_numpy_matrix(G) # Note: the order is G.nodes().
    adjacency_name = 'adjacency_' + noise + '.txt'
    np.savetxt(adjacency_name, A, fmt = '%1i', delimiter = ' ')

    ####################
    # Threshold matrix #
    ####################
    l = []
    for v in G.nodes():
        d_v = G.degree[v]
        l.append(np.random.randint(low = 1, high = d_v + 2, size = (c)))

    T = np.asarray(l)
    threshold_name = 'threshold_' + noise + '.txt'
    np.savetxt(threshold_name, T, fmt = '%1.3f', delimiter = ' ')

    #################
    # Weight matrix #
    ################# 
    W = np.random.rand(c,c)
    W = W/W.sum(axis=1)[:,None]
    weight_name = 'weight_' + noise + '.txt'
    np.savetxt(weight_name, W, fmt = '%1.10f', delimiter = ' ')

    ########################
    # Configuration matrix #
    ########################
    B = np.random.choice([0, 1], size = (n, c), p=[0.5, 0.5])
    # B = np.random.randint(low = 0, high = 2, size = (n, c)) # Note: exclusion on the high value
    configuration_name = 'configuration_' + noise + '.txt'
    np.savetxt(configuration_name, B, fmt = '%1.3f', delimiter = ' ')

    ##################
    # Density vector #
    ##################
    d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value
    density_name = 'density_' + noise + '.txt'
    np.savetxt(density_name, d, fmt = '%1.3f', delimiter = ' ')


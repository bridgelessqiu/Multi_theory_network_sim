import mms.threshold as mt
from scipy import sparse
import sys
import time
import mms.utility as mu
import networkx as nx
import numpy as np
import sys

"""
A and B are sparse
mt.correlate_threshold_density(A, B, T, d, k)
mt.correlate_threshold_count(A, B, T, W, k, r = 0)
mt.isolate_threshold_count(A, B, T, k, r = 0)

"""

if __name__ == "__main__":
    # command line args
    network_type = sys.argv[1] # real, synthetic
    exp_type = sys.argv[2] # increase_edges, increase_nodes
    network_name = sys.argv[3] # bio, bio2, ...
    number_of_contagions = sys.argv[4]

    # the number of contagions
    c = int(number_of_contagions)

    # the number of iterations
    num_iter = 200

    # parameters for synthetic networks
    low = 1000 # the size of the smallest network
    high = 50000 # the size of the largerst network
    step = 3000
    
    # percentage of infection
    p_infection = 0.5

    # list of p
    p_list = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.21, 0.23, 0.25]

    # list of all functions
    list_of_func = [mt.correlate_threshold_density, mt.correlate_threshold_count, mt.isolate_threshold_count]
    f = list_of_func[1] # The function to select
    
    if network_type == 'real':
        if exp_type == 'increase_nodes':
            # directoy
            dirc = 'datasets/'+ network_type + '/' + exp_type + '/' + network_name + '/'

            # file path
            f_path = network_name + ".edges"

            # read in the graph
            G = nx.read_edgelist(dirc + f_path, data=False)

            # graph info
            n = G.number_of_nodes()
            m = G.number_of_edges()
            print("n = {}, m = {}".format(n, m))

            ####################
            # Adjacency matrix #
            ####################
            A = nx.to_numpy_matrix(G) # Note: the order is G.nodes().

            ####################
            # Threshold matrix #
            ####################
            l = []
            for v in G.nodes():
                d_v = G.degree[v]
                l.append(np.random.randint(low = 1, high = d_v + 2, size = (c)))

            T = np.asarray(l)

            #################
            # Weight matrix #
            ################# 
            W = np.random.rand(c,c)
            W = W/W.sum(axis=1)[:,None]

            ########################
            # Configuration matrix #
            ########################
            B = np.random.choice([0, 1], size = (n, c), p=[1 - p_infection, p_infection])

            ##################
            # Density vector #
            ##################
            d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value

            # The starting processing time
            start = time.process_time()

            # Funciton call
            f(A, B, T, W, k = num_iter)

            # The ending processing time
            end = time.process_time()

            print("The time in seconds for {} is {}".format(network_name, end-start))
            print()
        
        elif exp_type == 'increase_edges':
            noise_levels = ['0', '0.01', '0.03', '0.05', '0.07', '0.09', '0.11', '0.13', '0.15', '0.17', '0.19', '0.21', '0.23', '0.25']

            # directoy
            dirc = 'datasets/'+ network_type + '/' + exp_type + '/' + network_name + '/'

            for noise in noise_levels:
                # file path
                f_path = network_name + '_' + noise + '_g2.edges'

                # read in the graph
                G = nx.read_edgelist(dirc + f_path, data=False)

                # the number of vertices
                n = G.number_of_nodes()

                ####################
                # Adjacency matrix #
                ####################
                A = nx.to_numpy_matrix(G) # Note: the order is G.nodes()

                ####################
                # Threshold matrix #
                ####################
                l = []
                for v in G.nodes():
                    d_v = G.degree[v]
                    l.append(np.random.randint(low = 1, high = d_v + 2, size = (c)))

                T = np.asarray(l)

                #################
                # Weight matrix #
                ################# 
                W = np.random.rand(c,c)
                W = W/W.sum(axis=1)[:,None]

                ########################
                # Configuration matrix #
                ########################
                B = np.random.choice([0, 1], size = (n, c), p=[1 - p_infection, p_infection])

                ##################
                # Density vector #
                ##################
                d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value

                start = time.process_time() # Start the time count

                # Function call
                f(A, B, T, W, k = num_iter)
                end = time.process_time() # End the time count

                print("The time in seconds for the network:{} and the percentage: {} is {}".format(network_name, noise, end-start))
                print()

    elif network_type == 'synthetic':
        if exp_type == 'increase_nodes':
            # directoy
            dirc = 'datasets/'+ network_type + '/' + exp_type + '/'

            for n in range(low, high, step):
                f_path = 'er_' + str(n) + '.edges'

                # read in the graph
                G = nx.read_edgelist(dirc + f_path, data=False)

                # Network info
                n = G.number_of_nodes()
                m = G.number_of_edges()
                print("n = {}, m = {}\n".format(n, m))

                ####################
                # Adjacency matrix #
                ####################
                A = nx.to_numpy_matrix(G) # Note: the order is G.nodes()
                A = sparse.csr_matrix(A)

                ####################
                # Threshold matrix #
                ####################
                l = []
                for v in G.nodes():
                    d_v = G.degree[v]
                    l.append(np.random.randint(low = 1, high = d_v + 2, size = (c)))
                T = np.asarray(l)

                #################
                # Weight matrix #
                #################
                W = np.random.rand(c,c)
                W = W/W.sum(axis=1)[:,None]

                ########################
                # Configuration matrix #
                ########################
                B = np.random.choice([0, 1], size = (n, c), p=[1.0 - p_infection, p_infection])
                B = sparse.csr_matrix(B)

                ##################
                # Density vector #
                ##################
                d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value

                start = time.time() # The starting time

                # The function call
                f(A, B, T, W, k = num_iter)

                end = time.time() # The ending time

                print("RE network of size {} takes {} seconds".format(n, end - start))
                print()
        elif exp_type == 'increase_edges':
            # directoy
            dirc = 'datasets/'+ network_type + '/' + exp_type + '/'

            for p in p_list:
                # The path
                f_path = 'er_' + str(p) + '.edges'

                # read in the graph
                G = nx.read_edgelist(dirc + f_path, data=False)

                # Network info
                n = G.number_of_nodes()
                m = G.number_of_edges()
                print("n = {}, m = {}\n".format(n, m))

                ####################
                # Adjacency matrix #
                ####################
                A = nx.to_numpy_matrix(G, dtype=int) # Note: the order is G.nodes()

                ####################
                # Threshold matrix #
                ####################
                l = []
                for v in G.nodes():
                    d_v = G.degree[v]
                    l.append(np.random.randint(low = 1, high = d_v + 2, size = (c)))
                T = np.asarray(l)

                #################
                # Weight matrix #
                #################
                W = np.random.rand(c,c)
                W = W/W.sum(axis=1)[:,None]

                ########################
                # Configuration matrix #
                ########################
                B = np.random.choice([0, 1], size = (n, c), p=[1.0 - p_infection, p_infection])

                ##################
                # Density vector #
                ##################
                d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value

                start = time.process_time() # Start the time count

                # Funciton call
                f(A, B, T, W, k = num_iter)

                end = time.process_time() # End the time count

                print("The ER graph with p = {} takes {} seconds".format(p, end-start))
                print()


                








            












import mms.threshold as mt
from guppy import hpy
from scipy import sparse
import sys
import time
import numpy as np
import sys

if __name__ == "__main__":
    # Experimental type
    type = 'edge'
    
    # Network size
    n = 200000

    # Average degree
    ave_deg = 100

    # Number of contagions
    c = 3

    # Number of iterations
    num_iter = 150

    # percentage of infection
    infection_dens = 0.3
    
    # functions
    list_of_func = [mt.correlate_threshold_density, mt.correlate_threshold_weight, mt.isolate_threshold_count]
    f = list_of_func[1]
    
    # Netowrk path name:
    path ="large/" + type + "/ER_" + str(n) + '_' + str(ave_deg) + ".npz"
    
    # The adjacency matrix
    A = sparse.load_npz(path)

    # Number of edges
    m = A.count_nonzero()/2

    # Network info
    print("The network has {} vertices and {} edges".format(n, m))

    # The degree array 
    deg = A.getnnz(1)

    # Threshold matrix
    l = []
    for i in range(n):
        d_i = deg[i]
        l.append(np.random.randint(low = 1, high = d_i + 2, size = (c)))
    T = np.asarray(l)

    # Weight matrix
    W = np.random.rand(c,c)
    W = W/W.sum(axis=1)[:,None]

    # Configuraiton matrix
    #B = np.random.choice([0, 1], size = (n, c), p=[1 - infection_dens, infection_dens])
    #B = sparse.csr_matrix(B)
    B = sparse.random(n, c, density = infection_dens, format='csr', data_rvs=np.ones)

    # Density vector
    d = np.random.randint(low = 1, high = c+1, size = (c, 1)) # Note: exclusion on the high value

    # Starting time
    start = time.time()

    # Funciton call
    f(A, B, T, W, k = num_iter)

    # The ending processing time
    end = time.time()

    print("RE network of size {} takes {} seconds".format(n, end - start))
    print("------------------------------\n")
    
    # Keep track of memory usage
    h = hpy()
    print(h.heap())

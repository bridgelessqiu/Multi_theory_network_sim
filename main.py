import mms.threshold as mt
import mms.utility as mu
import numpy as np
import sys

"""
correlate_threshold_density(A, B, T, d, k)

correlate_threshold_count(A, B, T, W, k, r = 0)

correlate_threshold_fraction(A, B, T, k, r = 0)

isolate_threshold_count(A, B, T, k, r = 0)
"""

if __name__ == "__main__":
    # command line args
    network_type = sys.argv[1]
    exp_type = sys.argv[2]
    network_name = sys.argv[3]
    
    # directoy
    dirc = 'datasets/'+ network_type + '/' + exp_type + '/' + network_name
    
    # Load all the data
    A = np.loadtxt(dirc+'/adjacency.txt')
    B = np.loadtxt(dirc+'/configuration.txt')
    T = np.loadtxt(dirc+'/threshold.txt')
    W = np.loadtxt(dirc+'/weight.txt')
    d = np.loadtxt(dirc+'/density.txt', ndmin=2)

    mt.correlate_threshold_count(A, B, T, W, k = 10)
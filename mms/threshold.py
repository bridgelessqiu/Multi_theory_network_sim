"""
Author: Zirou Qiu
Last modfied: 09/16/2020 
Description: 
    This module consists of simulations of the spread of multiple 
    contigions on a single network under the threshold model.
"""

#--------------------------- Imports ------------------------------#
# import networkx as nx
import numpy as np
import mms.utility as mu


#----------------------- Funciton Defintions ----------------------#

def isolate_threshold_count(A, B, T, k, r = 0, hub = False, clust = False):
    """
        Description
        -----------
        This function simulate the spread of multiple contigions on a single network
        where each contagion has 2 states (0 or 1).
        Contagions are not interrelated.

        Parameters
        ----------
        A: numpy matrix, int {0, 1}
          The adjacency matrix of G.

        B: numpy matrix, int {0, 1}
          The initial configuration matrix where $B_{vj}$ is the state value of 
          vertex v for contagion j.

        T: numpy matrix, int
          The threshold matrix where $T_{vj}$ is the threshold of vertex v for
          contagion j.

        k: int
          The number of system iterations
        
        r: float, optional
          The recovery probability. In each iteration, each vertex has a probability 
          r changing the state to 0 for each contigion.

        hub: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          hubs contract the contagions
        
        clust: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          vertices with the highest clustering coefficients contract the contagions.

        Returns
        -------
        B: numpy matrix
          The final configuration
    """

    # Check if A is indeed an adjacency matrix
    if not (mu.is_symmetric(A)):
        raise ValueError("The adjacency matrix is not symmetric")
    
    # Make all 1s along the diagonal of A (since we are considering the closed neighborhood)
    np.fill_diagonal(A, 1)
    
    # The recovery probability
    recovery = False
    if r != 0:
        recovery = True

    # The main loop
    for i in range(k):
        # matrix operation
        B_last = B
        B = A*B_last - T

        # update states
        B[B >= 0] = 1
        B[B < 0] = 0

        # If a recovery probability is set
        if recovery:
            B[np.random.rand(*B.shape) < r] = 0

        # if fixed point
        if np.array_equal(B, B_last):
            print("A fixed point is reached at iteration {}".format(i))
            return B

    return B

###########################################################################################

def correlate_threshold_count(A, B, T, W, k, r = 0, hub = False, clust = False):
    """
        Description
        -----------
        This function simulate the spread of multiple contigions on a single network
        where each contagion has 2 states (0 or 1).
        Contagions are not interrelated.

        Parameters
        ----------
        A: numpy matrix, int {0, 1}
          The adjacency matrix of G.

        B: numpy matrix, int {0, 1}
          The initial configuration matrix where $B_{vj}$ is the state value of 
          vertex v for contagion j.

        T: numpy matrix, int
          The threshold matrix where $T_{vj}$ is the threshold of vertex v for
          contagion j.

        W: numpy matrix, float [0, 1]
          The weight matrix where $W_{ij}$ is the weight of contagion j w.r.t 
          contagion i

        k: int
          The number of system iterations
        
        r: float, optional
          The recovery probability. In each iteration, each vertex has a probability 
          r changing the state to 0 for each contigion.

        hub: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          hubs contract the contagions
        
        clust: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          vertices with the highest clustering coefficients contract the contagions.

        Returns
        -------
        B: numpy matrix
          The final configuration
    """

    # Check if A is indeed an adjacency matrix
    if not (mu.is_symmetric(A)):
        raise ValueError("The adjacency matrix is not symmetric")
    
    # Make all 1s along the diagonal of A (since we are considering the closed neighborhood)
    np.fill_diagonal(A, 1)
    
    # The recovery probability
    recovery = False
    if r != 0:
        recovery = True

    # Take the transpose of the weight matix
    W = np.transpose(W)

    # The main loop
    for i in range(k):
        # matrix operation
        B_last = B
        B = A * B_last * W - T

        # update states
        B[B >= 0] = 1
        B[B < 0] = 0

        # If a recovery probability is set
        if recovery:
            B[np.random.rand(*B.shape) < r] = 0

        # if fixed point
        if np.array_equal(B, B_last):
            print("A fixed point is reached at iteration {}".format(i))
            return B

    return B

###########################################################################################

def correlate_threshold_fraction(A, B, T, k, r = 0, hub = False, clust = False):
    """
        Description
        -----------
        This function simulate the spread of multiple contigions on a single network
        where each contagion has 2 states (0 or 1).
        Contagions are not interrelated.

        Parameters
        ----------
        A: numpy matrix, int {0, 1}
          The adjacency matrix of G.

        B: numpy matrix, int {0, 1}
          The initial configuration matrix where $B_{vj}$ is the state value of 
          vertex v for contagion j.

        T: numpy matrix, float [0, 1]
          The threshold matrix where $T_{vj}$ is the threshold (floats) of vertex v for
          contagion j.

        k: int
          The number of system iterations
        
        r: float, optional
          The recovery probability. In each iteration, each vertex has a probability 
          r changing the state to 0 for each contigion.

        hub: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          hubs contract the contagions
        
        clust: boolean, optional
          If set true, prints out the indices of the iterations for which the top 10 
          vertices with the highest clustering coefficients contract the contagions.

        Returns
        -------
        B: numpy matrix
          The final configuration
    """

    # Check if A is indeed an adjacency matrix
    if not (mu.is_symmetric(A)):
        raise ValueError("The adjacency matrix is not symmetric")
    
    # Make all 1s along the diagonal of A (since we are considering the closed neighborhood)
    np.fill_diagonal(A, 1)

    # Number of vertices
    n = np.shape(A)[0]

    # The inverse of diagonal degree matrix
    d = np.reciprocal(
        (A * np.ones((n, 1), dtype = 'float')).A1
    )
    D = np.diag(d)

    # The recovery probability
    recovery = False
    if r != 0:
        recovery = True

    # The main loop
    for i in range(k):
        # matrix operation
        B_last = B
        B = D * A * B_last - T

        # update states
        B[B >= 0] = 1.0
        B[B < 0] = 0.0

        # If a recovery probability is set
        if recovery:
            B[np.random.rand(*B.shape) < r] = 0.0

        # if fixed point
        if np.array_equal(B, B_last):
            print("A fixed point is reached at iteration {}".format(i))
            return B

    return B

###########################################################################################

def correlate_threshold_density(A, B, T, d, k):
    """
        Description
        -----------
        This function simulate the spread of multiple contigions on a single network
        where each contagion has 2 states (0 or 1).
        Contagions are not interrated as described by the second model.

        Parameters
        ----------
        A: numpy matrix, int {0, 1}
          The adjacency matrix of G.

        B: numpy matrix, int {0, 1}
          The initial configuration matrix where $B_{vj}$ is the state value of 
          vertex v for contagion j.

        T: numpy matrix, float [0, 1]
          The threshold matrix where $T_{vj}$ is the threshold (floats) of vertex v for
          contagion j.

        d: numpy matrix, int
          The density vector where $d_j$ is the density of the contagion $c_j$

        k: int
          The number of system iterations

        Returns
        -------
        B: numpy matrix
          The final configuration
    """
    # Compute the reciprocal 
    d_bar = np.transpose( np.reciprocal(d.astype(float)) )

    # The number of contagions
    c =  np.shape(T)[1]

    # k * 1 ones
    one = np.ones((c, 1), dtype = 'float')

    for i in range(k):
      B_last = B

      # Compute the totoal number of infections of each vertex
      b = B * one

      # Compute M
      M = b * d_bar
      M[M >= 1.0] = 1.0
      M[M < 1.0] = 0.0
      
      B = A * M - T

      # update states
      B[B >= 0.0] = 1.0
      B[B < 0.0] = 0.0

      # if fixed point
      if np.array_equal(B, B_last):
          print("A fixed point is reached at iteration {}".format(i))
          return B
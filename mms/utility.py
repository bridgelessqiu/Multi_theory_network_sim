"""
Author: Zirou Qiu
Last modfied: 09/16/2020 
Description: 
    This module consists of utility funcitons
"""

import numpy as np

# Check if a matrix is symmetric
def is_symmetric(A, rtol=1e-05, atol=1e-08):
    # Source: https://stackoverflow.com/a/42913743
    return np.allclose(A, A.T, rtol=rtol, atol=atol)
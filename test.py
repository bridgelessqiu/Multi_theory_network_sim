import mms.threshold as mt
import mms.utility as mu
import numpy as np

if __name__ == "__main__":
    A = np.matrix([[1, 2],[3, 4]])
    C = np.matrix([[4, 4],[4, 4]])
    B = np.matrix([[2, 2], [2, 2]])
    B[A<C] = 0
    print(B)
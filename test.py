import mms.threshold as mt
import mms.utility as mu
import numpy as np

def test(x):
    return x + 1 

def test2(fun, i):
    return fun(i)

if __name__ == "__main__":
    print(test2(test, 5))
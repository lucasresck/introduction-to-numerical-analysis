import numpy as np

def fixed_point(g, a, b, x_0, tol, k):
    err = tol + 1
    x = x_0
    while err > tol:
        x_old = x
        x = g(x)
        err = k/(1-k)*np.abs(x-x_old)
    return x
    
def g(x):
    return np.log(15 - np.log(x))
    
result = fixed_point(g, 1, 3, 1, 10**-5, 0.1)
print(result, g(result))


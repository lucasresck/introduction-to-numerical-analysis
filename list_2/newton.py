import numpy as np

def newton(h, dh, x_0, tol):
    err = tol + 1
    x = x_0
    while err > tol:
        x_old = x
        x = x - h(x)/dh(x)
        err = np.abs(x-x_old)/np.abs(x_old)
    return x

def h(x):
    return np.log(15 - np.log(x)) - x
    
def dh(x):
    return 1/(x*np.log(x) - 15*x) - 1
    
result = newton(h, dh, 2, 10**-8)
print(result, h(result))

import numpy as np

def collocation(k):
    z = [[0]*k]
    for i in range(1,k+1):
        z[0][i-1]= (i-(1/2))/k
    z = np.array(z)
    return z

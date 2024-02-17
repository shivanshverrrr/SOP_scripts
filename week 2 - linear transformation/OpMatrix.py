import Haar
import numpy as np
import math

alpha = 1
bound = 1

def felem(k):
    #returns elements of matrix F
    val = ((k+1)**(alpha+1)) - (2*(k**(alpha+1))) + ((k-1)**(alpha+1))
    return val

def fmatrix(n):
    #forms the matrix F
    f = np.zeros((n,n))
    for i in range(0,n):
        for j in range(0,n):
            if j-i>0:
                f[i,j]= felem(j-i)
            elif j==i:
                f[i,j]= 1
            f[i,j] = f[i,j]/(math.gamma(alpha+2)*(n**alpha))
    return f

def opmatrix(n,hmat):
    #bound**alpha = 1 here
    return (bound**alpha)*np.matmul(hmat,np.matmul(fmatrix(n),np.linalg.inv(hmat)))
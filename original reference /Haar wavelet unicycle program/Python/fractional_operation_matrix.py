import numpy as np
from element_of_matrix_F import element_of_matrix_F 
from math import gamma
from HaarMatrix import HaarMatrix
from casadi import *

def fractional_operation_matrix(n,alpha,bound,H):
    F = np.zeros((n,n))

    for i in range(1,n+1):
        for j in range(1,n+1):
            if j-i>0:
                F[i-1][j-1] = element_of_matrix_F(j-i,alpha)

            elif j==i:
                F[i-1][j-1] = 1
            else:
                F[i-1][j-1]=0
            
            F[i-1][j-1] = F[i-1][j-1]/(gamma(alpha+2)*(n**alpha))

    palpha =  (bound ** alpha) * np.matmul(H,np.matmul(F,np.linalg.inv(H)))
    
    return palpha
# print("original = ")
# print(fractional_operation_matrix(16,1,1,HaarMatrix(16))[:,[15]])

# print(fractional_operation_matrix(16,1,1,HaarMatrix(16)) @ HaarMatrix(16))
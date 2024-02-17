import Haar
import OpMatrix
import numpy as np
import math
import matplotlib.pyplot as plt
from casadi import *
opti = casadi.Opti()

#start and stop positions
x1,y1 = 0,0
x2,y2 = 1,1

#this can be used for linear transformations as required.
x1_1,y1_1 = 0,0.5
x2_2,y2_2 = 0.5,1

k = 32

hmat = Haar.hmatrix(k)
opmat = OpMatrix.opmatrix(k,hmat)

#these functions are used to transform the numpy matrices to casadi supported matrics
def hcol(k,t): #transforms the numpy row to a column matrix in casadi
    A = MX.zeros(k,1)
    a = Haar.hcolumn(t,k)
    for i in range(k):
        A[i] = a[i]
    return A
def palpha(k):
    A = MX.zeros(k,k)
    a = OpMatrix.opmatrix(k,hmat)
    for i in range(k):
        for j in range(k):
            A[i,j] = a[i,j]
    return A
def conv_hmat(k):
    A = MX.zeros(k,k)
    for i in range(0,k):
        for j in range(0,k):
            A[i,j] = hmat[i,j]
    return A

#Initializing the variables we need to solve for
C1 = opti.variable(1,k)
C2 = opti.variable(1,k)
C3 = opti.variable(1,k)  
C4 = opti.variable(1,k)

#derivatives represent velocity or the control variables
def xdot(t):
    return C1@hcol(k,t) # @ means matrix multiplication
def ydot(t):
    return C2@hcol(k,t)
def vel(t):
    return C3@hcol(k,t)
def theta(t):
    return C4@hcol(k,t)

#system variables/position
def x(t):
    return (C1@palpha(k))@hcol(k,t) + x1
def y(t):
    return (C2@palpha(k))@hcol(k,t) + y1

#The galerkian equations. Formed for every i and t (no. of i and t are equal)
def gal1(i,t):
    return Haar.haarfun(i,t) * ( xdot(t) - ( vel(t)@cos(theta(t)) ) )
def gal2(i,t):
    return Haar.haarfun(i,t) * ( ydot(t) - ( vel(t)@sin(theta(t)) ) )

A = MX.zeros(1,1) 
B = MX.zeros(1,1)
coloc = Haar.colloc(k)

G1 = MX.zeros(k,1) #column vectors
G2 = MX.zeros(k,1)

for i in range(0,k):
    A[0,0] = 0
    B[0,0] = 0
    for j in range(0,k):
        A[0,0] = A[0,0] + gal1(i,coloc[j])
        B[0,0] = B[0,0] + gal2(i,coloc[j])
    G1[i,0] = (1/k)*(A[0,0])
    G2[i,0] = (1/k)*(B[0,0])

#defining the cost function
A = MX.zeros(1,k)
B = MX.zeros(1,1)
cv_hmat = conv_hmat(k)
for i in range(k):
    A[0,i] = C3@cv_hmat[:,i]
    B[0,0] = B[0,0] + A[0,i]**2
cost = (1/k)*B

#Equality constraints
XEQ = ((C1@palpha(k))@cv_hmat[:,k-1]) + x1 - x2
YEQ = ((C2@palpha(k))@cv_hmat[:,k-1]) + y1 - y2

#velocity and angle constraints
x = (C1@palpha(k))@cv_hmat + x1
y = (C2@palpha(k))@cv_hmat + y1
v = C3@cv_hmat
theta = C4@cv_hmat

# opti.minimize(cost)
# opti.subject_to(G1 == 0)
# opti.subject_to(G2 == 0)
# opti.subject_to(XEQ == 0)
# opti.subject_to(YEQ == 0)

# #obstacle constraints
# #class circle for quicker initializaton of obstacles
# class circle():
#     def __init__(self,cx,cy,r):
#         self.cx = cx
#         self.cy = cy
#         self.r = r
#     def eq(self):
#         return (x-self.cx)**2 + (y-self.cy)**2 - self.r**2
#     def patch(self):
#         return plt.Circle((self.cx,self.cy),self.r,fill = False)

# #O B S T A C L E S
# ob1 = circle(0.5,0.5,0.1)
# opti.subject_to(ob1.eq()>0)

# #velocity constraints
# opti.subject_to(v<1.5)
# opti.subject_to(v>-1.5)

# #angle constraints
# opti.subject_to(theta < 3.14)
# opti.subject_to(theta > 0)

# opti.solver('ipopt')
# sol = opti.solve()

# c1 = sol.value(C1)
# c2 = sol.value(C2)
# c3 = sol.value(C3)
# c4 = sol.value(C4)

# #extracting values of x and y
# x = np.zeros((1,k))
# y = np.zeros((1,k))
# for i in range(k):
#     x[0,i] = x1_1 + (x2_2-x1_1)*np.matmul(np.matmul(c1,opmat),hmat[:,i]) 
#     y[0,i] = y1_1 + (y2_2-y1_1)*np.matmul(np.matmul(c2,opmat),hmat[:,i])

# #-------------------------------------------------------------------------#
# print(x)
# print(y)
# b = sol.value(B)
# print("cost = " + str(b/k))

# #drawing graphs
# figure,axes = plt.subplots()
# plt.plot(x,y,marker = 'x',color = 'black')

# axes.add_patch(ob1.patch())

# plt.plot()
# plt.show()

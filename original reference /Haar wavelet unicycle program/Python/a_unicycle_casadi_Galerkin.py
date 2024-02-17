import sys 
sys.path.append('../Helpers') # import the helpers folder

from casadi import *
from HaarMatrix import HaarMatrix
from fractional_operation_matrix import fractional_operation_matrix
import math
import numpy as np
from collocation import collocation
import matplotlib.pyplot as plt

opti = Opti()

k = 32
alpha = 1
b = 1
x1_0 = 0
x2_0 = 0
x1_1=1
x2_1=1

H = HaarMatrix(k)
palpha = fractional_operation_matrix(k, alpha, b, H)

# get the obstacles
obstacles = [
    {
        "x": 0.4,
        "y": 0.5,
        "r": sqrt(0.1)
    },
    {
        "x": 0.8,
        "y": 1.5,
        "r": sqrt(0.1)
    }
]

no_of_obstacles = len(obstacles)

x = opti.variable(1, 4*k)  # 4*k unknowns

C1T = x[0, 0:k]
C2T = x[0, k:2*k]
C3T = x[0, 2*k:3*k]
C4T = x[0, 3*k:4*k]



aabc = np.zeros((1, k))
for i in range(1, k+1):
    aabc[0, i-1] = 2**(-1*math.floor(np.log(i)/np.log(2)))
bbcc = horzcat(1, aabc)


first_DE = MX.zeros(1, k)
second_DE = MX.zeros(1, k)


uu = MX.zeros(k, k)
vv = MX.zeros(k, k)
ww = MX.zeros(k, k)

A = MX.zeros(k, 1)
B = MX.zeros(k, 1)
C = MX.zeros(k, 1)

# calculate the constraints using galerkins method
for i in range(0, k):
    for j in range(0, k):
        uu[i, j] = H[i,j] * cos(C3T @ H[:, [j]])
        vv[i, j] = H[i, j] * sin(C3T @ H[:, [j]])
        A[i] = A[i] + uu[i,j]
        B[i] = B[i] + vv[i,j]
    first_DE[0, i] = bbcc[0, i] *( (C1T @ palpha[:, [k-1]]) -(1/k)*vv[i,j]*A[i, 0])
    second_DE[0,i] = bbcc[0,i]* ( (C2T @ palpha[:, [k-1]])  -(1/k)*vv[i,j]*B[i,0])


F1 = first_DE
F2 = second_DE

# define x y in terms of variables
xxx=(C1T @ palpha  @ H)  + x1_0
yyy=(C2T @ palpha  @ H)  + x2_0

# define the obstacle equations
obstacle_eqn = []
for i in range(no_of_obstacles):
    x_center_circle = obstacles[i]['x']
    y_center_circle = obstacles[i]['y']
    r_circle = obstacles[i]['r']
    obstacle_eqn.append(-((xxx-x_center_circle)**2 + (yyy-y_center_circle)**2 - r_circle**2))


F11=(C1T @ palpha @ H[:,[k-1]]  + x1_0)  -   x1_1   #   equality constraint i.e. boundary condtion for x1 at 0  
F12=(C2T @ palpha @ H[:,[k-1]]  + x2_0)   -  x2_1   #   equality constraint i.e. boundary condtion for x1 at 0  

# define the cost function
xxa = 0
for i in range(0,k):
    aaa = C3T @ H[:,[i]]
    bbb = C4T @ H[:,[i]]
    xxa += aaa**2 


Cost = xxa/k

# use casadi to solve the problem
opti.minimize(Cost)
opti.subject_to(F1==0)
opti.subject_to(F2==0)


opti.subject_to(F11==0)
opti.subject_to(F12==0)


for i in range(no_of_obstacles):
    opti.subject_to(obstacle_eqn[i]<=0)

# try solving the problem if the solution is feasible then plot it on the graph with the obstacles
try:
    opti.solver('ipopt')
    sol = opti.solve()
    sol.value(x)

    c1t=sol.value(C1T)
    c2t=sol.value(C2T)
    c3t=sol.value(C3T)
    c4t=sol.value(C4T)


    tt = collocation(k)
    xxb = np.zeros((1,k))
    xxc = np.zeros((1,k))
    for i in range(0,k):
        xxb[0,i]= c1t@palpha@H[:,[i]]  +  x1_0
        xxc[0,i]= c2t@palpha@H[:,[i]]  +  x2_0

    figures,axes = plt.subplots()
    for i in range(no_of_obstacles):
        circle = plt.Circle((obstacles[i]['x'],obstacles[i]['y']),obstacles[i]['r'],fill = False)
        axes.add_artist(circle)
    axes.set_aspect(1)


    plt.plot(xxb[0],xxc[0],'-x')
    plt.xlim(0,1.5)
    plt.ylim(0,2.5)
    plt.xlabel("x axis")
    plt.ylabel("y axis")

    plt.show()
except KeyboardInterrupt:
    print('Exiting.')
except: # else exit the program and just plot the obstacles
    print("Unable to solve the problem")
    figures,axes = plt.subplots()
    for i in range(no_of_obstacles):
        circle = plt.Circle((obstacles[i]['x'],obstacles[i]['y']),obstacles[i]['r'],fill = False)
        axes.add_artist(circle)
    axes.set_aspect(1)
    plt.xlim(0,1.5)
    plt.ylim(0,2.5)
    plt.xlabel("x axis")
    plt.ylabel("y axis")

    plt.show()


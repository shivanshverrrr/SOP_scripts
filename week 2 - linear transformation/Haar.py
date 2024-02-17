#contains methods related to operations on haar functions and matrices.
import numpy as np
from matplotlib import pyplot as plt
import math

def haarfun(n,t):
    #returns value of a haar function given any value of n and t.

    #defining parent functions to calculate further functions.
    def mother(t):
        if 0<=t and t<0.5:
            return 1
        elif 0.5<=t and t<1:
            return -1
        else:
            return 0
    def scaling(t):
        if 0<=t and t<1:
            return 1
        else:
            return 0

    if n == 0:
        return scaling(t)
    if n == 1:
        return mother(t)

    if math.floor(math.log2(n+1)) == math.ceil(math.log2(n+1)):
        j = math.floor(math.log2(n+1)) - 1
    else:
        j = math.floor(math.log2(n+1))
    k = n - 2**j

    return mother((2**j)*t - k)

def colloc(k): 
    #returns a list of points collocated around 0.5
    a = []
    for i in range(1,k+1): #includes 1 and excludes k+1
        a.append((i-0.5)/k)
    a = np.array(a)
    return a

def graph_hf(n):
    #use as debugging tool. graphs a given haar function for a value of n
    #works well with small sizes of n, else increase num.
    x = np.linspace(start = 0,stop = 1,num = 500, endpoint = False) #exclued endpoints.
    #num is how many samples to generate
    y = []
    for i in x:
        y.append(haarfun(n,i))
    y = np.array(y)
    #print(y)
    plt.plot(x,y)
    plt.show()

def hcolumn(t,k): 
    #NOTE: even though the function reads hcolumn, the return type is a
    #row vector. It's just easier to implement most methods using numpy
    #when the return elements is a row vector. the function is named hcolumn
    #to stay true to the notation used in the papers.

    #returns haar column vector for a given value of t with k elements.
    #the ith element corresponds to the (i-1)th haar function
    #the values of t are calculated using collocation around 0.5.
    a = []
    for i in range(0,k):
        a.append(haarfun(i,t))
    a = np.array(a)
    return a

def hmatrix(k):
    #reutrns a square matrix where each columns corresponds to the value of the first k
    #haar functions at a given time t. The value of the time is obtained using 
    #collocation. 
    a = np.zeros((k,k))
    tlist = colloc(k)
    col = 0
    for i in tlist:
        a[:,col] = hcolumn(i,k) 
        col = col + 1
    return a
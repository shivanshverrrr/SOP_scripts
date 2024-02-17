import numpy as np

def elements_of_Haar_natrix(k,j,n):
    t = (2*j+1)/(2*n)
    p =  np.floor(np.log(k)/np.log(2))
    up = 2**p
    q = k-up+1

    low = (q-1)/up
    mid = (q-0.5)/up 
    high = q/up 

    if (t>=low) and  (t<mid):
        value = 2**(p/2)
    elif t>=mid and t<high:
        value = -2**(p/2)
    else:
        value = 0
    
    return value
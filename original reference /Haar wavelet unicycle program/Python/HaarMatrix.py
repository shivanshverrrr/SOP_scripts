import numpy as np

# Discrete wavelets and pertubation theory - http://lib.ysu.am/articles_art/3fa03dd41bdd763fca71fce5487d88e2.pdf
# https://en.wikipedia.org/wiki/Haar_wavelet

def HaarMatrix(N):
    n = np.floor(np.log(N)/np.log(2))
    if 2**n!=N:
        print("Error: size "+str(N)+" is not a multiple of power 2")
        return -1

    if N>2:
        h = HaarMatrix(N/2) # recursively find the previous haar matrix
    else:
        return np.array([[1,1],[1,-1]]) # base case
    
    h_upper = np.kron(h,[1,1]) # calculate the upper part of the Haar Matrix

    h_lower = np.kron(np.identity(int(N/2)),[1,-1]) # calculate the lower part of the Haar Matrix

    h = np.vstack((h_upper,h_lower)) # combine both the parts

    return h


# print(HaarMatrix(8))
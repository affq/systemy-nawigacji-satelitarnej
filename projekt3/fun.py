import math
import numpy as np

a = 6378137
e2 = 0.00669438002290

def Np(B):
    N = a/(1-e2*(np.sin(B)**2))**0.5
    return N

def Hirvonen(X,Y,Z):
    r = (X**2 + Y**2)**0.5
    B = math.atan(Z/(r*(1-e2)))
    
    while 1:
        N = Np(B)
        H = r/np.cos(B) - N
        Bst = B
        B = math.atan(Z/(r*(1-(e2*(N/(N+H))))))    
        if abs(Bst-B)<(0.00001/206265):
            break
    L = math.atan2(Y,X)
    N = Np(B)
    H = r/np.cos(B) - N
    return B, L, H

def Rneu(phi, lamb):
    R = np.array([[-np.sin(phi)*np.cos(lamb), -np.sin(lamb), np.cos(phi)*np.cos(lamb)],
                    [-np.sin(phi)*np.sin(lamb), np.cos(lamb), np.cos(phi)*np.sin(lamb)],
                    [np.cos(phi), 0, np.sin(phi)]])
    return R
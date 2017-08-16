from __future__ import division
import sys
#import matplotlib.pyplot as plt
import numpy as np
import math as m
import random as rd
homeModules = 'C:\\Users\\Claire\\Documents\\GitHub\\HumanEvol'
sys.path.append(homeModules)

def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = np.multiply(maxVal, np.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2)))
    return v
    
def ReprodBiAllelic(freq):
    p = freq[0,:]
    q = freq[1,:]
    pPrime = p**2 + np.multiply(p,q)
    qPrime = q ** 2 + np.multiply(p,q)
    freqPrime = np.vstack((pPrime,qPrime))
    freqPrime = freqPrime/np.sum(freqPrime)
    return freqPrime
    
def ReproTriAllelic(freq):
    p1 = freq[0,:]
    p2 = freq[1,:]
    p3 = freq[2,:]
    p1Prime = p1**2 + np.multiply(p1,p2) + np.multiply(p1,p3)
    p2Prime = p2 ** 2 + np.multiply(p1,p2) + np.multiply(p2,p3)
    p3Prime = p3**2 + np.multiply(p1,p3) + np.multiply(p2,p3)
    freqPrime = np.vstack((p1Prime,p2Prime,p3Prime))
    freqPrime = freqPrime/np.sum(freqPrime)
    return freqPrime    
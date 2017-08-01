# -*- coding: utf-8 -*-
"""
Claire Guerin - 31/07/2017
"""

# run in command line: reset -f

from __future__ import division
import matplotlib.pyplot as plt
import math
import numpy as np
import time


def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = maxVal * math.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2))
    return v
    
nGen = 10000
nStrategies = 100
m = 0.1 # upmigration capacity

fmax = 0.2    

p = np.empty([nGen, nStrategies, nStrategies]) # the population is monomorphic
p[1,:,:] = 1 # the population is monomorphic at t=0
f1 = 0.1 # proportion of individuals in class 1 at t=0

#xr = rep(NA,nGen) # strategy of the resident population
#xr[1] = 1 # strategy of the resident population at t=0
xstrat = np.linspace(0, 1, nStrategies) # strategies

# class 1 viability
mu1 = 0.5
sigma1 = 0.1

# class 2 viability
mu2 = 0.3
sigma2 = 0.2
max2 = 0.8

startTime = time.clock()

for j in range(nStrategies):
        
    xr = xstrat[j]
    
    for k in range(nStrategies):
            
        xm = xstrat[k]

        for i in range(nGen - 1):
            
            if (fmax <= f1):
                
                m21 = 0 # proportion of C2 that upmigrate
                m12 = f1 - fmax # proportion of C1 that downmigrate
                
            else:
                        
                m21 = m * (fmax - f1) / (1 - f1)
                m12 = 0

            pr = p[i, j, k]
            rs1 = pr * f1 * (1 - m12) * Viability(xr, mu1, sigma1) + pr * (1 - f1) * m21 * Viability(xr, mu2, sigma2, max2)
            rs2 = pr * (1 - f1) * (1 - m21) * Viability(xr, mu2, sigma2, max2) + pr * f1 * m12 * Viability(xr, mu1, sigma1)
        
            ms1 = (1 - pr) * f1 * (1 - m12) * Viability(xm, mu1, sigma1) + (1 - pr) * (1 - f1) * m21 * Viability(xm, mu2, sigma2, max2)
            ms2 = (1 - pr) * (1 - f1) * (1 - m21) * Viability(xm, mu2, sigma2, max2) + (1 - pr) * f1 * m12 * Viability(xm, mu1, sigma1)
            
            divisor = rs1 + ms1 + rs2 + ms2
            
            f1 = (rs1 + ms1) / divisor

            p[i + 1, j, k] = (rs1 / divisor) ** 2 + (rs2 / divisor) ** 2 + ms1 / divisor * rs1 / divisor + ms2 / divisor * rs2 / divisor
    
endTime = time.clock()-startTime

#nCheckPoints = 3
#mat = np.empty([nStrategies, nStrategies, nCheckPoints])
#
#plt.figure(1)
#
#for n in range(nCheckPoints):
#    
#    gen = [1000,10000,100000][n]
#    mat[:,:,gen] = p[:,:,gen]
#    plt.subplot(nCheckPoints,1,n)
#    plt.contourf(mat[:,:,n]) 
#    
#stDev = np.std(mat,3)
#nonEquilibrium = stDev > 0
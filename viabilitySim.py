# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# run in command line: reset -f

import math
import numpy as np
import time

def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = maxVal * math.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2))
    return v
    
nGen = [1000, 10000, 100000]
nRun = len(nGen)
runDict = {}
nStrategies = 100
m = 0.1 # upmigration capacity
f1 = 0.1 # proportion of individuals in class 1 at t=0
fmax = 0.2

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

for l in range(nRun):
    
    runTime = nGen[l]
    p = np.empty([runTime, nStrategies, nStrategies]) # the population is monomorphic
    p[1,:,:] = 1 # the population is monomorphic at t=0
    
    for j in range(nStrategies):
        
        xr = xstrat[j]
    
        for k in range(nStrategies):
            
            xm = xstrat[k]

            for i in range(runTime - 1):
            
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

    runDict[runTime] = p
    del p
    
endTime = time.clock()-startTime
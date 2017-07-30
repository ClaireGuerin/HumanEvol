# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 12:07:34 2017

@author: Claire
"""

import matplotlib.pyplot as plt
import numpy as np

runDict = {1000:[],10000:[],100000:[]}
allRuns = sorted(runDict.keys())
nRuns = len(allRuns)
nStrategies = 100

mat = np.empty([nStrategies, nStrategies, nRuns])

plt.figure(1)

for n in range(nRuns):
    
    key = allRuns[n]
    mat[:,:,n] = runDict[key][-1]
    plt.subplot(3,1,n)
    plt.contourf() 
    
stDev = np.std(mat,3)
nonEquilibrium = stDev > 0
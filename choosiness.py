# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:49:04 2017

@author: Claire
"""

# run in command line: reset -f

from __future__ import division
import matplotlib.pyplot as plt
import math
import numpy as np
import time


def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = maxVal * np.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2))
    return v

def ReproductiveSuccess(choosiness):
    f = - choosiness^2 + choosiness + 3/4 
    return(f)
    
r = 1/10 # recombination rate
t = 100 # number of generations

ht = np.empty([t,4]) # haplotype frequencies over time (4 haplotypes)
at = np.empty([t,2]) # resident alleles frequencies over time (p and q)
dt = np.empty([t,1]) # linkage disequilibrium D between X and Phi over time

ht[0,0] = 1/10
ht[0,1] = 6/10
ht[0,2] = 1/10
ht[0,3] = 2/10

at[0,0] = ht[0,0] + ht[0,1]
at[0,1] = ht[0,0] + ht[0,2]

dt[0,0] = ht[0,0] * ht[0,3] - ht[0,1] * ht[0,2]

for gen in range(1,t):
    ht[gen,0] = (1 - r) * ht[gen-1,0] + r * at[gen-1,0] * at[gen-1,1]
    ht[gen,1] = (1 - r) * ht[gen-1,1] + r * at[gen-1,0] * (1 - at[gen-1,1])
    ht[gen,2] = (1 - r) * ht[gen-1,2] + r * (1 - at[gen-1,0]) * at[gen-1,1]
    ht[gen,3] = (1 - r) * ht[gen-1,3] + r * (1 - at[gen-1,0]) * (1 - at[gen-1,1])
    
    at[gen,0] = ht[gen,0] + ht[gen,1]
    at[gen,1] = ht[gen,0] + ht[gen,2]

    dt[gen,0] = (1 - r) * dt[gen-1,0]

plt.plot(ht) 
plt.plot(at)
plt.plot(dt)
    
nStrategies = 100
m = 0.1 # upmigration capacity

fmax = 0.2    

p = np.empty([nGen, nStrategies, nStrategies]) # the population is monomorphic
p[0,:,:] = 1 # the population is monomorphic at t=0
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

class Population():
    
    def __init__(self):
        self.haplotypeFrequencies
        self.recombinationRate
        self.mutationRate
        self.mutationStep
        
        self.freqClass1
        
        self.p = self.haplotypeFrequencies[0] + self.haplotypeFrequencies[1]
        self.q = self.haplotypeFrequencies[0] + self.haplotypeFrequencies[2]

        self.dietTrait = np.empty(2)
        self.choosinessTrait = np.empty(2)
        
        self.muDiet = np.empty(2)
        self.sigmaDiet = np.empty(2)
        self.maxDiet = np.empty(2)
        
        self.linkageDisequilibrium = self.haplotypeFrequencies[0] * self.haplotypeFrequencies[3] - self.haplotypeFrequencies[1] * self.haplotypeFrequencies[2]



    def nutritionSelection(self,xmut,xres):
        survivalClass1 = self.freqClass1 * np.array[self.p,1-self.p] * Viability(self.dietTrait, self.muDiet[0], self.sigmaDiet[0], self.maxDiet[0])
        survivalClass2 = (1 - self.freqClass1) * np.array[self.p,1-self.p] * Viability(self.dietTrait, self.muDiet[1], self.sigmaDiet[1], self.maxDiet[1])
        
        
        
        
        
        
        
        residentClass1 =  + self.p * (1 - self.freqClass1) * self.migration21 * Viability(self.residentDiet, muClass2, sigmaClass2, maxClass2)
        residentClass2 =  + pr * f1 * m12 * Viability(xr, mu1, sigma1)
        
        mutantClass1 = (1 - pr) * f1 * (1 - m12) * Viability(xm, mu1, sigma1) + (1 - pr) * (1 - f1) * m21 * Viability(xm, mu2, sigma2, max2)
        mutantClass2 = (1 - pr) * (1 - f1) * (1 - m21) * Viability(xm, mu2, sigma2, max2) + (1 - pr) * f1 * m12 * Viability(xm, mu1, sigma1)
            
        meanViability = rs1 + ms1 + rs2 + ms2
            
        f1 = (rs1 + ms1) / meanViability
    
        p[i + 1, j, k] = (rs1 / meanViability) ** 2 + (rs2 / meanViability) ** 2 + ms1 / meanViability * rs1 / meanViability + ms2 / meanViability * rs2 / meanViability
    

for j in range(nStrategies):
        
    xr = xstrat[j]
    phir = xstrat[j]
    
    for k in range(nStrategies):
            
        xm = xstrat[k]
        phim = xstrat[k]

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
            
            meanViability = rs1 + ms1 + rs2 + ms2
            
            f1 = (rs1 + ms1) / meanViability

            p[i + 1, j, k] = (rs1 / meanViability) ** 2 + (rs2 / meanViability) ** 2 + ms1 / meanViability * rs1 / meanViability + ms2 / meanViability * rs2 / meanViability
    
endTime = time.clock()-startTime

nCheckPoints = 3
mat = np.empty([nStrategies, nStrategies, nCheckPoints])

plt.figure(1)

for n in range(nCheckPoints):
    
    gen = np.linspace(0,nGen,nCheckPoints)[n]
    mat[:,:,gen] = p[:,:,gen]
    plt.subplot(nCheckPoints,1,n)
    plt.contourf(mat[:,:,n]) 
    
stDev = np.std(mat,3)
nonEquilibrium = stDev > 0





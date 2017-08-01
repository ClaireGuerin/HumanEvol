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
    """ g1 = f(BC), g2 = f(Bc), g3 = f(bC), g4 = f(bc) 
    where B is adaptation to diet and C is choosiness"""
    
    def __init__(self):
        self.haplotypes = []
        self.recombinationRate
        self.mutationRate
        self.mutationStep
        
        self.freqClass1
        self.classes = [self.freqClass1,1-self.freqClass1]
        
        self.p = self.haplotypeFrequencies[0] + self.haplotypeFrequencies[1]
        self.q = self.haplotypeFrequencies[0] + self.haplotypeFrequencies[2]

        self.dietTrait = []
        self.choosinessTrait = []
        
        self.muDiet = [] # mu class 1, mu class 2
        self.sigmaDiet = [] # sigma class 1, sigma class 2
        self.maxDiet = [] # max survival class 1, max survival class 2
        
        self.linkageDisequilibrium = self.haplotypeFrequencies[0] * self.haplotypeFrequencies[3] - self.haplotypeFrequencies[1] * self.haplotypeFrequencies[2]



    def nutritionSelection(self):
        
        haplotypes = np.transpose(np.array(self.haplotypes).reshape(1,4))
        classes = np.array(self.classes).reshape(1,2)
        
        typeFrequencies = haplotypes * classes
        
        x = np.tile(np.repeat([self.dietTrait,1-self.dietTrait],2).reshape(1,4),2)
        mu = np.tile(self.muDiet,4).reshape(4,2)
        sigma = np.tile(self.sigmaDiet,4).reshape(4,2)
        maxi = np.tile(self.maxDiet,4).reshape(4,2)
        
        survival = typeFrequencies * Viability(x,mu,sigma,maxi)
        meanViability = np.sum(survival)
            
        newFreqClass1 = np.sum(survival,0)[0] / meanViability

        residentFreq = np.sum(survival[[0,1],:],0) / meanViability
        mutantFreq = np.sum(survival[[2,3],:],0) / meanViability
    
        #newPResident = np.sum(residentFreq ** 2 + residentFreq * mutantFreq) # this is for reroduction i.e. frequencies at the next generation
        
        return(newPResident, newFreqClass1)
    
        
        
    def socialMigration(self, fmax):
        
        f1 = self.freqClass1
        
        if (fmax <= f1):
            m21 = 0 # proportion of C2 that upmigrate
            m12 = f1 - fmax # proportion of C1 that downmigrate
                
        else:
            m21 = m * (fmax - f1) / (1 - f1)
            m12 = 0
            
        haplotypes = np.transpose(np.array(self.haplotypes).reshape(1,4))
       
        return [m21,m12]
        
        
        
        
        

for j in range(nStrategies):
        
    xr = xstrat[j]
    phir = xstrat[j]
    
    for k in range(nStrategies):
            
        xm = xstrat[k]
        phim = xstrat[k]

        for i in range(nGen - 1):
            
            

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





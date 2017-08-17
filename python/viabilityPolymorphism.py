# -*- coding: utf-8 -*-
"""
Claire Guerin - 31/07/2017
"""

# run in command line: reset -f

from __future__ import division
import sys
import numpy as np
import itertools as it
import time
import scipy.io as sio
homeModules = 'C:\\Users\\Claire\\Documents\\GitHub\\HumanEvol\\python'
sys.path.append(homeModules)
import UtilitiesViability as ut

dirname = 'C:\\Users\\Claire\\Dropbox\\MEME\\Montpellier\\AdaptiveDynamicsStratification\\viabilityRun'

# Import data with polymorphism values
#pData = sio.loadmat('%s\pGen0-2000.mat' % (dirname))['prob']
#startingMat = pData[-1,:,:]

# PARAMETERS 
NbOfGenerations = 1000
NbOfStrategies = 100  
NbOfAlleles = 2

# viability class param
Mus = [0.5,0.3] # up / down
Sigmas = [0.1,0.2] # up / down
Maxis = [1.0,0.8] # up / down

migUp = 0.2  
pUpTown = 0.1  
       

xstrat = np.linspace(0, 1, NbOfStrategies) # strategies
allStrat = np.tile(xstrat,NbOfAlleles).reshape(NbOfAlleles,len(xstrat))
loopOn = list(it.product(*allStrat))
NbOfCombinations = len(loopOn)
    
p = np.empty([NbOfGenerations,NbOfAlleles,NbOfCombinations])
startMono = 0.98
p[0,:,:] = np.tile(np.vstack((startMono,np.repeat(1-startMono,NbOfAlleles-1))),NbOfCombinations)

startTime = time.clock()

for alleleComb in range(NbOfCombinations):

    
    Xs = np.array(loopOn[alleleComb])
    probaInit = np.reshape(p[0,:,alleleComb],[NbOfAlleles,1])
    freqMat = np.multiply(probaInit,[pUpTown,1-pUpTown])
    
    for gen in range(1,NbOfGenerations):
        
        # 1: VIABILITY SELECTION
        x = np.repeat(Xs,2).reshape(NbOfAlleles,2)
        mu = np.tile(Mus,NbOfAlleles).reshape(NbOfAlleles,2)
        sigma = np.tile(Sigmas,NbOfAlleles).reshape(NbOfAlleles,2)
        maxi = np.tile(Maxis,NbOfAlleles).reshape(NbOfAlleles,2)
        
        survival = ut.Viability(x,mu,sigma,maxi)  
        survival = np.multiply(freqMat,survival)
        survival = survival / np.sum(survival)        
                
        # 2: MIGRATION
        migration = survival
        migration[:,0] = migration[:,0] + migUp * migration[:,1]
        migration[:,1] = (1 - migUp) * migration[:,1]
        migration = migration/np.sum(migration)
        
        # 3: REPRODUCTION 
        
        if NbOfAlleles > 2:
            reproduction = ut.ReproTriAllelic(migration)
        else:
            reproduction = ut.ReprodBiAllelic(migration)
        
        # 4: REGULATION
        
        reproduction[:,0] = reproduction[:,0] - survival[:,1] * migUp
        reproduction[:,1] = reproduction[:,1] + survival[:,1] * migUp
        reproduction = reproduction/np.sum(reproduction)
        freqMat = reproduction
        
        pPrime = np.sum(reproduction,1)
        p[gen,:,alleleComb] = pPrime
    
endTime = time.clock()-startTime

myDict = {'prob': p, 'strat': loopOn}
sio.savemat('%s\\viabEvol.mat' % dirname, myDict)

#import operator
#getX = operator.itemgetter(0)
#xMut = list(map(getX, loopOn))
#xMut > 0.5
#xMut < 0.5

#nBits = 100-62
#bits = np.linspace(lastGen,200000,nBits + 1).astype(int)
#
#for i in range(len(bits)):
#    
#    cutoff = range(bits[i],bits[i+1])
#    mat2save = p[cutoff,:,:]
#    myDict = {'prob': mat2save}
#    sio.savemat('%s\pGen%s-%s.mat' % (dirname, str(bits[i]), str(bits[i+1])),myDict)

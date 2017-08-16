from __future__ import division
from multiprocessing import Pool
from contextlib import closing
import sys
#import time
import itertools as it
#import matplotlib.pyplot as plt
import numpy as np
import random as rd
import scipy.io as sio
homeModules = '/n/home06/cguerin/src/humanevol'
sys.path.append(homeModules)
from Parameters import *
import Utilities as ut
from Population import *


dirname = '/n/regal/debivort_lab/claire/humanevol'

t = 100 # number of generations

#p = np.empty([t, NbOfStrategies, NbOfStrategies])
#q = np.empty([t, NbOfStrategies, NbOfStrategies])
p = 1
q = 1 # the population is monomorphic at t=0

residentDiets = np.linspace(0,1,NbOfStrategies)
residentPhis = np.linspace(0,1,NbOfStrategies)

recombinationRate = [.5,.75,1] # no recomb, some recomb, full recomb

fertilityBenef = [0.01,0.1]   

heritD = [.25,.5,.75] # heritability of trait D
heritA = [.25,.5,.75] # heritability of trait A


m0 = [.1,.5]         
h0 = [.1,.5]
seeding = list(range(NbOfStrategies))

allParams = [seeding,recombinationRate,list(residentDiets),list(residentPhis),heritA,heritD,m0,h0,fertilityBenef]
allParamCombinations = list(it.product(*allParams))

def runEvolution(parcomb):
    
    parcomb = int(parcomb)
    seedingInst = allParams[parcomb,0]
    recombRInst = allParams[parcomb,1]
    resDietInst = allParams[parcomb,2]
    resPhiInst = allParams[parcomb,3]
    heritAInst = allParams[parcomb,4]
    heritDInst = allParams[parcomb,5]
    m0Inst = allParams[parcomb,6]
    h0Inst = allParams[parcomb,7]
    bInst = allParams[parcomb,8]

    rd.seed(seedingInst)
    
    evolveMat = np.empty([4,16,t])
    
    popEvolve = Population()
    
    popEvolve.r = recombRInst
    popEvolve.dietResidentMutant = [resDietInst,resDietInst]
    popEvolve.phiResidentMutant = [resPhiInst,resPhiInst]
    popEvolve.Hattractive = heritAInst
    popEvolve.Hdominant = heritDInst
    popEvolve.socialMigration[0] = m0Inst
    popEvolve.hypergyny[0] = h0Inst
    popEvolve.b = bInst

    evolveMat[:,:,0] = popEvolve.startMatrix
    popEvolve.nutritionSelection(popEvolve.startMatrix)
    
    for gen in range(1,t):
        # MIGRATION
        popEvolve.socialMigration()
        
        # MUTATION: does mutant replace resident? does a mutation occur?
        ut.updateTrait(popEvolve.dietResidentMutant)
        ut.updateTrait(popEvolve.phiResidentMutant)
        
        # REPRODUCTION
        upperclass = popEvolve.migMat[:,range(8)]
        lowerclass = popEvolve.migMat[:,range(8,16)]
        RepUpperClass = pUpTown*popEvolve.matingGenotype(upperclass)
        RepLowerClass = (1-pUpTown)*popEvolve.matingGenotype(lowerclass)
        allFrequencies = np.hstack((RepUpperClass,RepLowerClass))
        popEvolve.repMat = allFrequencies/np.sum(allFrequencies)
        
        # REGULATION
        popEvolve.downRegulation()
        evolveMat[:,:,gen] = popEvolve.regMat
        
        # SELECTION
        popEvolve.nutritionSelection(popEvolve.regMat)
    
    myDict = {'evolGenPhen': evolveMat}
    sio.savemat('%s/runseed%sr%sdiet%sphi%s-ha%shd%smbas%shbas%sb%s.mat' % (dirname, str(seedingInst), str(recombRInst), str(resDietInst), str(resPhiInst),str(heritAInst),str(heritDInst), str(m0Inst),str(h0Inst),str(bInst)),myDict)
        

if __name__ == '__main__':
    # start 4 worker processes
    with closing(Pool(processes=32)) as pool:

        pool.map(runEvolution, allParamCombinations)      
        pool.terminate
        
        
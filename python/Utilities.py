from __future__ import division
import sys
import os
#import matplotlib.pyplot as plt
import numpy as np
import math as m
import random as rd
homeModules = 'C:\\Users\\Claire\\Documents\\GitHub\\HumanEvol'
sys.path.append(homeModules)
from Parameters import *

def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = np.multiply(maxVal, np.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2)))
    return v

def ReproductiveSuccess(choosiness):
    f = - choosiness^2 + choosiness + 3/4 
    return f
 
def hapTransform(haplotypes):
    formattedHaplotypes = np.transpose(np.array(haplotypes).reshape(1,4))
    return(formattedHaplotypes)
 
def freqTransform(frequency):
    if isinstance(frequency,int):
        
        eltSize = 1
        
    elif isinstance(frequency,list):
        
        eltSize = len(frequency)

    formattedFrequencies = np.array(frequency).reshape(1,eltSize)
    return formattedFrequencies
    
    
def genMutation(geneticTrait):
    # genetic trait is either dietResMut or phiResMut
    
    if rd.uniform(0,1) < 0.5:
        mutationDirection = geneticTrait[0] + stepMutation
  
    else:
        mutationDirection = geneticTrait[0] - stepMutation

    if mutationDirection < 0:
        mutationDirection =  m.fabs(mutationDirection) 
    
    elif mutationDirection > 1:
        mutationDirection = 1 - 1 - mutationDirection  
                    
    geneticTrait[1] = mutationDirection
    return geneticTrait
    
    
    
def updateTrait(trait):
    
    if trait[1] > 0.999:
        trait = [trait[1],0]
    else:
        trait = trait
        
    if rd.uniform(0,1.0) < mutationRate and trait[1] == 0:
        # there is a mutation
        trait = genMutation(trait)
    else:
        trait = trait
        
    return trait
    
    
def testExtinction(matrix):
    if not np.count_nonzero(matrix):
        matrix = matrix
    else:
        matrix = matrix/np.sum(matrix)
    return matrix
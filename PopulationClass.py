# -*- coding: utf-8 -*-
"""
Created on Fri Jul 07 14:05:43 2017

@author: Claire
"""

from __future__ import division
import itertools as it
import numpy as np
import matplotlib.pyplot as plt

mutationRate = 0.0005
mutationStep = 0.01

mutation = 0 # no mutation at this step; if mutation: 1

class Population():
    
    def __init__(self):
        self._typeNames = ["Attractive", "DietQuality", "Choosy", "Dominant", "Female", "UpTown"]
        self._typeFreqs = [0, 0, 0, 0, 0.5, 0.1] 
        # frequency of each type, including sex ratio & probability to be from the higher stratum
        
        self._nTypes = len(self._typeNames)
        self._allComb = list(it.product([0, 1], repeat = self._nTypes))
        
        
        self._deltaMig = 0.1
        
        self._pResidentUpMigration = np.array([[0.7,0.3],[0.2,0.8]]) # transition matrix: chances to migrate from one social class to the other for the resident dominance level
        self._pMutantUpMigration = self._pResidentUpMigration + np.tile([-self._deltaMig, self._deltaMig], (2, 1)) # transition matrix: chances to migrate from one social class to the other for the mutant dominance level
        
        self._dietRequirement = np.array([0.75, 0.75]) # diet quality requirement 
        
    def Combinations(self):
        
        self._nComb = len(self._allComb)
        self._allProb = np.empty([self._nComb, len(self._allComb[0])])
        
        for i in range(self._nComb):
            for j in range(self._nTypes):
                if self._allComb[i][j]:
                    self._allProb[i][j] = self._typeFreqs[j]
                else:
                    self._allProb[i][j] = 1 - self._typeFreqs[j]

        return np.sum(self._allProb, axis = 1)

    
    def Migration(self):
        
        self._pD_class = [self._typeFreqs[2] * (1 - self._pUpTown) * (1 - self._pFemale), self._typeFreqs[2] * self._pUpTown * (1 - self._pFemale)]  # proportion of males of resident dominance in each social class
        self._pD_classMig = np.dot(self._pD_class, self._pUpMigration) # proportion of males of resident dominance in each social class after migration
        
    def Selection(self):
        self._meanDiet = np.array([0.75, 0.25])
        self._stdDiet = np.array([0.25, 0.25]) 
        # mean & std diet requirement of viability selection in habitats "High social stratum" & "Low social stratum" resp.
        self._pSurvival = np.exp(-(self._dietRequirement - self._meanDiet) ** 2 / (2 * self._stdDiet ** 2)) # survival probability according to diet requirement in each habitat
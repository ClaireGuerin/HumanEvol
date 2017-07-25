# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 17:19:34 2017

@author: Claire
"""
from __future__ import division


import scipy.stats as scist
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
#import math as m

# Linkage disequilibrium male choosiness - local adaptation to nutrition
# Local adaptation and male choosiness are genetically determined by two bi-allelic loci B & C, with alleles B & b and C & c.

# Female attractivity and male social dominance are two phenotypes alpha and delta with two morphs a0 & a1 and d0 & d1 respectively. They are genetically and environmentally determined, with an heritability of h_A and h_D respectively.

# if h_A = h_D = 1, we can write the system as a system of 4 loci A, B, C and D.
# Let the loci be independant. Thus, r = 1/2 for every loci combination. 


haplotypes = ['ABCD','ABCd','ABcD','ABcd','AbCD','AbCd','AbcD','Abcd','aBCD','aBCd','aBcD','aBcd','abCD','abCd','abcD','abcd']
               
def cross(haplo):
    return[x for x in it.combinations_with_replacement(haplo,2)]
    
crossings = cross(haplotypes)
len(crossings)

               
#crossings = it.combinations(haplotypes,2)

# Haplotypes frequencies B & C

r = 1/10

t = 100 #generations

ht = np.empty([t,4]) # haplotype frequencies over time (8 haplotypes)
at = np.empty([t,2]) # dominant alleles frequencies over time (B and C)
dt = np.empty([t,1]) # linkage disequilibrium D over time

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
    
    
#f1_prime = f1_f * f1_m + 1/2 * (f1_f * (f2_m + f3_m) + f1_m * (f2_f + f3_f)) + 1/2 * (1 - r) * (f1_f * f4_m + f1_m * f4_f) + 1/2 * r * (f2_f * f3_m + f2_m * f3_f)
#
#f2_prime = f2_f * f2_m + 1/2 * (f2_f * (f1_m + f4_m) + f2_m * (f1_f + f4_f)) + 1/2 * (1 - r) * (f2_f * f3_m + f2_m + f3_f) + 1/2 * r * (f1_f * f4_m + f1_m * f4_f)
#
#f3_prime = f3_f * f3_m + 1/2 * (f3_f * (f1_m + f4_m) + f3_m * (f1_f + f4_f)) + 1/2 * (1 - r) * (f2_f * f3_m + f2_m + f3_f) + 1/2 * r * (f1_f * f4_m + f1_m * f4_f)
#
#f4_prime = f4_f * f4_m + 1/2 * (f4_f * (f2_m + f3_m) + f4_m * (f2_f + f3_f)) +  1/2 * (1 - r) * (f1_f * f4_m + f1_m * f4_f) + 1/2 * r * (f2_f * f3_m + f2_m * f3_f)

# Phenotypes 0,0 ; 0,1; 1,0; 1,1

hA = 3/5
hD = 1/5

mutationRate = 0.0005
degreeOfMutation = 0.1

class IndividualProfile():
    """Class defining an individual according to: 
        - its phenotype (attractivity false/true, dominance false/true)
        - its genotype (mutant diet quality requirement false/true, mutant choosiness false/true)
        - its sex (female false/true)
        - its social class (high false/true)"""
        
    def __init__(self, attractivity, dominance, choosiness, diet, sex, social_class):
        self.attractivity = attractivity
        self.dominance = dominance
        self.choosiness = choosiness
        self.diet = diet
        self.sex = sex
        self.social_class = social_class
        
    def PhenProd(self):
        #with which probability can it produce each type of phenotype ?
        
        if self.attractivity & self.dominance:
            self.a0d0 = hA * hD
            self.a0d1 = hA * (1 - hD)
            self.a1d0 = (1 - hA) * hD
            self.a1d1 = (1 - hA) * (1 - hD)
            
        if self.attractivity & self.dominance == False:
            self.a0d0 = hA * (1 - hD)
            self.a0d1 = hA * hD
            self.a1d0 = (1 - hA) * (1 - hD)
            self.a1d1 = (1 - hA) * hD

        if self.attractivity == False & self.dominance:
            self.a0d0 = (1 - hA) * hD
            self.a0d1 = (1 - hA) * (1 - hD)
            self.a1d0 = hA * hD
            self.a1d1 = hA * (1 - hD)
            
        if self.attractivity == False & self.dominance == False:
            self.a0d0 = (1 - hA) * (1 - hD)
            self.a0d1 = (1 - hA) * hD
            self.a1d0 = hA * (1 - hD)
            self.a1d1 = hA * hD
            
#    def Fertility(self):
#        
#        if self.sex:

##########################################################################
##########################################################################

mutationRate = 0.0005
mutationStep = 0.01

mutation = 0 # no mutation at this step; if mutation: 1

class Population():
    
    def __init__(self):
        self._typeNames = ["Attractive","DietQuality","Choosy","Dominant","Female","UpTown"]
        self._typeFreqs = [0,0,0,0,0.5,0.1] 
        # frequency of each type, including probability to be female & probability to be from the higher stratum
        
        self._allComb = list(it.product([0, 1], repeat=len(self._typeNames)))
        
        self._deltaMig = 0.1
        
        self._pResidentUpMigration = np.array([[0.7,0.3],[0.2,0.8]]) # transition matrix: chances to migrate from one social class to the other for the resident dominance level
        self._pMutantUpMigration = self._pResidentUpMigration + np.tile([-self._deltaMig, self._deltaMig], (2, 1)) # transition matrix: chances to migrate from one social class to the other for the mutant dominance level
        
        self._dietRequirement = np.array([0.75, 0.75]) # diet quality requirement  
    
    def Migration(self):
        
        self._pD_class = [self._typeFreqs[2] * (1 - self._pUpTown) * (1 - self._pFemale), self._typeFreqs[2] * self._pUpTown * (1 - self._pFemale)]  # proportion of males of resident dominance in each social class
        self._pD_classMig = np.dot(self._pD_class, self._pUpMigration) # proportion of males of resident dominance in each social class after migration
        
    def Selection(self):
        self._meanDiet = np.array([0.75, 0.25])
        self._stdDiet = np.array([0.25, 0.25]) 
        # mean & std diet requirement of viability selection in habitats "High social stratum" & "Low social stratum" resp.
        self._pSurvival = np.exp(-(self._dietRequirement - self._meanDiet) ** 2 / (2 * self._stdDiet ** 2)) # survival probability according to diet requirement in each habitat
        

x = np.linspace(0,1,100)
meanX = 0.75
stdX = 0.25
surv = np.exp(-(x - meanX) ** 2 / (2 * stdX ** 2))

plt.plot(x,surv)

allClasses = list(it.product([0, 1], repeat=6))



pHypergyny = 0.01 # chances to be considered for inter-class mating for the resident female attractivity level --> offspring will thus be raised in higher social stratum
pFecundity = 0.5 # chances to successfully reproduce after mating for the resident female attractivity level

hA = 3/5 # probability for parent attractivity to be transmitted to offspring
hD = 1/5 # probability for parent dominance to be transmitted to offspring



offspring = np.empty([16,16])
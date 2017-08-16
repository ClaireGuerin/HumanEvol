# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:49:04 2017

@author: Claire
"""

# run in command line: reset -f

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import time
import Utilities as ut
    
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

p = np.empty([t, nStrategies, nStrategies]) # the population is monomorphic
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
        
        self.r # recombination rate
        self.mutationRate
        self.mutationStep
        
        self.pUpTown = .1 # starting frequency of uptown population
        self.sexRatio = .5
        self.probaB = 1
        self.probaC = 1
        self.probad0 = .5
        self.probaa0 = .5

        self.dietResidentMutant = []
        self.phiResidentMutant = []
        
        self.muDiet = [] # mu class 1, mu class 2
        self.sigmaDiet = [] # sigma class 1, sigma class 2
        self.maxDiet = [] # max survival class 1, max survival class 2
        
        self.s = 0
        self.e = 0
        
        self.socialMigration = [] # m0, m, m'
        self.hypergyny = [] # h0, h, h'

        dummyCol = np.repeat(.25,8)*self.sexRatio
        colMat = np.hstack((dummyCol*self.pUpTown,dummyCol*(1-self.pUpTown)))
        rowMat = np.transpose([1,0,0,0])
        self.startMatrix = np.multiply(rowMat,colMat) 
        # all genotypes & phenotypes/sex/social class
        
        self.LD = rowMat[0]*rowMat[3]-rowMat[1]*rowMat[2]



    def nutritionSelection(self,matrix):
        # matrix here should either be startMat or regMat
        
        x = np.repeat([self.dietResidentMutant[0],self.dietResidentMutant[1]],2).reshape(2,2)
        mu = np.tile(self.muDiet,2).reshape(2,2)
        sigma = np.tile(self.sigmaDiet,2).reshape(2,2)
        maxi = np.tile(self.maxDiet,2).reshape(2,2)
        
        survival = ut.Viability(x,mu,sigma,maxi)
        resSurv = np.vstack((np.resize(survival[0,0],[2,8]),np.resize(survival[0,1],[2,8])))
        mutSurv = np.vstack((np.resize(survival[1,0],[2,8]),np.resize(survival[1,1],[2,8])))
        
        survival = np.multiply(matrix,np.hstack((resSurv,mutSurv)))
        selMat = survival / np.sum(survival)
        
        return selMat
        



        
        
    def socialMigration(self, matrix):
        # matrix here should be selMat
        
        ms = self.socialMigration
        m0 = ms[0]
        m = ms[1]
        mPrime = ms[2]
        MUpd0 = m0
        MUpd1 = m0 * mPrime/m
        
        hs = self.hypergyny
        h0 = hs[0]
        h = hs[1]
        hPrime = hs[2]
        Ha0 = h0
        Ha1 = h0 * hPrime/h
        
        matrix[:,(0,2)] = matrix[:,(0,2)] + MUpd0 * matrix[:,(8,10)]
        matrix[:,(0,2)] = (1 - MUpd0) * matrix[:,(8,10)]
        matrix[:,(1,3)] = matrix[:,(1,3)] + MUpd1 * matrix[:,(9,11)]
        matrix[:,(1,3)] = (1 - MUpd1) * matrix[:,(9,11)]
        matrix[:,(4,5)] = matrix[:,(4,5)] + Ha0 * matrix[:,(12,13)]
        matrix[:,(4,5)] = (1 - Ha0) * matrix[:,(12,13)]
        matrix[:,(6,7)] = matrix[:,(6,7)] + Ha1 * matrix[:,(14,15)]
        matrix[:,(6,7)] = (1 - Ha1) * matrix[:,(14,15)]
        
        migMat = matrix/np.sum(matrix)

        return migMat
        
        
        
        
        
#        
#    def mating(self):
#        
#        haplotypes = self.freqMigration
#        r = self.recombinationRate
#        
#        highStratumMales = 1/2 * haplotypes[:,0]
#        lowStratumMales = 1/2 * haplotypes[:,1]
#
#        females = 1/2 * np.sum(haplotypes,1) # haplotype frequencies in the whole female pool
#
#        newHaplotypes = np.empty([4,4]) # haplotypes in 4 categories (classes and sex)
#        
#        newHaplotypes[0,:] = g1class
#
#        
#        haplotypes[0,:] = (1 - r) * ht[gen-1,0] + r * at[gen-1,0] * at[gen-1,1]
#        haplotypes[1,:] = (1 - r) * ht[gen-1,1] + r * at[gen-1,0] * (1 - at[gen-1,1])
#        haplotypes[2,:] = (1 - r) * ht[gen-1,2] + r * (1 - at[gen-1,0]) * at[gen-1,1]
#        haplotypes[3,:] = (1 - r) * ht[gen-1,3] + r * (1 - at[gen-1,0]) * (1 - at[gen-1,1])
#    
#        at[gen,0] = ht[gen,0] + ht[gen,1]
#        at[gen,1] = ht[gen,0] + ht[gen,2]
#
#        dt[gen,0] = (1 - r) * dt[gen-1,0]

    upclass = self.allprobs[:,range(8)]
    loclass = self.allprobs[:,range(8,16)]
                            
                            
                            
                            
                            
                            
                            
                            
    def available(self,matrix,PhiS):
        
        s = self.s
        e = self.e
        
        lm = self.latency[0]
        lf = self.latency[1]
        
        phip = PhiS[0]
        phi = PhiS[1]
        pa1 = np.sum(matrix[:,range(6,8)])/np.sum(matrix[:,range(4,8)]) 
        
        alpha_res = np.sqrt(phip) # to replace with mathematica solution for dr
        alpha_a1 = 1/(1 + (s*e*alpha_res*lf)/(1 - s*lf))
        alpha_a0 = 1/(1 + s*e*alpha_res*(1 - phip)*lf/(1 - s*lf))
        alpha_mut = 1/(1 + s*e*lm*(pa1*alpha_a1 + (1 - pa1)*alpha_a0*(1 - phi))/(1-s*lm))
        
        availabilities = np.array([alpha_res,alpha_mut,alpha_a0,alpha_a1])
        return availabilities
        
        
        
        
                            
        
    def phenoRep(self,matrix, availability,phi):
        #here, matrix is the submatrix for one class, one allele of choosiness with associated phi value, i.e. 2 rows of genotypes
        
        b = self.b
        
        a0 = availability[2]*(1-phi)
        a1 = availability[3]*(1+b) # b is the benefit obtained from mating with an attractive female
        
        phenoRep = np.vstack((np.multiply(matrix[:,0],matrix[:,4]*a0), 
                    np.multiply(matrix[:,0],matrix[:,5]*a0) + np.multiply(matrix[:,1],matrix[:,4]*a0),
                    np.multiply(matrix[:,0],matrix[:,6]*a1) + np.multiply(matrix[:,2],matrix[:,4]*a0),                 
                    np.multiply(matrix[:,0],matrix[:,7]*a1) + np.multiply(matrix[:,3],matrix[:,4]*a0),
                    np.multiply(matrix[:,1],matrix[:,5]*a0),
                    np.multiply(matrix[:,1],matrix[:,6]*a1) + np.multiply(matrix[:,2],matrix[:,5]*a0),
                    np.multiply(matrix[:,1],matrix[:,7]*a1) + np.multiply(matrix[:,3],matrix[:,5]*a0),
                    np.multiply(matrix[:,2],matrix[:,6]*a1),
                    np.multiply(matrix[:,2],matrix[:,7]*a1) + np.multiply(matrix[:,3],matrix[:,6]*a1),
                    np.multiply(matrix[:,3],matrix[:,7]*a1)))
        
        return phenoRep
                        
                        
        
        
        
        
                            
    def matingPhenotype(self,matrix):
        
        hi = self.Hattractive * self.Hdominant
        ha = self.Hattractive * (1-self.Hdominant)
        hd = (1-self.Hattractive) * self.Hdominant
        hdelta = (1-self.Hattractive) * (1-self.Hdominant) 
        
        heritp1 = np.array([hi,1/2*(ha+hi),1/2*(hd+hi),1/2*(hdelta+hi),ha,1/2*(hd+ha),1/2*(hdelta+ha),hd,1/2*(hdelta+hd),hdelta])
        heritp2 = np.array([ha,1/2*(hi+ha),1/2*(hdelta+ha),1/2*(hd+ha),hi,1/2*(hdelta+hi),1/2*(hd+hi),hdelta,1/2*(hd+hdelta),hd])
        heritp3 = np.array([hd,1/2*(hdelta+hd),1/2*(hi+hd),1/2*(ha+hd),hdelta,1/2*(hi+hdelta),1/2*(ha+hdelta),hi,1/2*(ha+hi),ha])
        heritp4 = np.array([hdelta,1/2*(hd+hdelta),1/2*(ha+hdelta),1/2*(hi+hdelta),hd,1/2*(ha+hd),1/2*(hi+hd),ha,1/2*(hi+ha),hi])
        
        herit = np.transpose(np.vstack((heritp1,heritp2,heritp3,heritp4)))
        
        phenRep = self.phenoRep(matrix)
        
        phenNew = np.empty([4,4])
        
        for phen in range(4):
            phenNew[phen,] = np.multiply(phenRep,herit[phen])
            
        return phenNew

        
        
        
        
        
        
    def matingGenotype(self,matrix,LD):
        #matrix here is the matrix for 1 social class only
        
        phip = self.phiResidentMutant[0]
        phi = self.phiResidentMutant[1]

        s = self.s
        e = self.e

        availability = self.available(self,phip,phi)
        phenMat = self.matingPhenotype(matrix)
        
        matFemales = .5*phenMat
        matMales = .5*phenMat
        
        sevec = np.repeat(s*e,4)
        alphavec = np.resize([availability[0],availability[1]],4)
        
        matMales = np.multiply(matMales,sevec)
        matMales = np.multiply(matMales,alphavec)
        
        g1prime = matMales[:,0]*matFemales[:,0]+.5*(matMales[:,0]*matFemales[:,1]+matMales[:,1]*matFemales[:,0])+.5*(matMales[:,0]*matFemales[:,2]+matMales[:,2]*matFemales[:,0])+.5*(1-r)*(matMales[:,0]*matFemales[:,3]+matMales[:,3]*matFemales[:,0])+.5*r*(matMales[:,1]*matFemales[:,2]+matMales[:,2]*matFemales[:,1])
        
        g2prime = np.multiply(matMales[:,1],matFemales[:,1])+.5*(np.multiply(matMales[:,0],matFemales[:,1])+np.multiply(matMales[:,1],matFemales[:,0]))+.5*(np.multiply(matMales[:,1],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,1]))+.5*(1-r)*(np.multiply(matMales[:,1],matFemales[:,2])+np.multiply(matMales[:,2],matFemales[:,1]))+.5*r*(np.multiply(matMales[:,0],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,0]))
        
        g3prime = np.multiply(matMales[:,2],matFemales[:,2])+.5*(np.multiply(matMales[:,0],matFemales[:,2])+np.multiply(matMales[:,2],matFemales[:,0]))+.5*(np.multiply(matMales[:,2],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,2]))+.5*(1-r)*(np.multiply(matMales[:,1],matFemales[:,2])+np.multiply(matMales[:,2],matFemales[:,1]))+.5*r*(np.multiply(matMales[:,0],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,0]))
        
        g4prime = np.multiply(matMales[:,3],matFemales[:,3])+.5*(np.multiply(matMales[:,2],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,2]))+.5*(np.multiply(matMales[:,1],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,1]))+.5*r*(np.multiply(matMales[:,1],matFemales[:,2])+np.multiply(matMales[:,2],matFemales[:,1]))+.5*(1-r)*(np.multiply(matMales[:,0],matFemales[:,3])+np.multiply(matMales[:,3],matFemales[:,0]))
        
        genNew = 0.5*np.vstack((g1prime,g2prime,g3prime,g4prime))
        matNew = np.hstack((genNew,genNew))
        matNew = matNew/np.sum(matNew)
        
        LDNew = (1-r)*LD
        
        return matNew, LDNew        
        
    
        
        
        
        
        
    def downRegulation(self, repmatrix, selmatrix):
        
        c2Sel = selmatrix[:,range(8,12)]
                          
        MupTot = self.socialMigration[0]*(1+self.socialMigration[2]/self.socialMigration[1])
        HypTot = self.hypergyny[0]*(1+self.hypergyny[2]/self.hypergyny[1])
        
        repmatrix[:,range(4)] = repmatrix[:,range(4)] - c2Sel * (MupTot + HypTot)
        repmatrix[:,range(8,12)] = repmatrix[:,range(8,12)] + c2Sel * (MupTot + HypTot)
               
        regMat = repmatrix/np.sum(repmatrix)
        
        return regMat
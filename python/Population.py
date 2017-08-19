# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:49:04 2017

@author: Claire
"""

from __future__ import division
import numpy as np
import random as rd
import Utilities as ut
from Parameters import * 

class Population():
    """ g1 = f(BC), g2 = f(Bc), g3 = f(bC), g4 = f(bc) 
    where B is adaptation to diet and C is choosiness"""
    
    def __init__(self):
        
        self.r = 0 # recombination rate
        
        self.s = survivalRate
        self.e = encounterRate
        self.latency = [latencyMale,latencyFemale]

        self.sexRatio = sexRatio
        self.pUpTown = pUpTown
        
        self.muDiet = muDiet 
        self.sigmaDiet = sigmaDiet
        self.maxDiet = maxDiet

        self.dietResidentMutant = []
        self.phiResidentMutant = []

        self.Hattractive = 0 # heritability trait A
        self.Hdominant = 0 # heritability trait D
        
        self.socMigration = [0,migrationParamLow,migrationParamHigh] # m0, m, m'
        self.hypergyny = [0,migrationParamLow,migrationParamHigh] # h0, h, h'
        
        self.b = 0 # fertility advantage for a1 females

        self.probaB = probaB
        self.probaC = probaC
        self.probad0 = probad0
        self.probaa0 = probaa0

        dummyCol = np.repeat(.25,8)*self.sexRatio
        colMat = np.hstack((dummyCol*self.pUpTown,dummyCol*(1-self.pUpTown)))
        rowMat = np.array([1,0,0,0]).reshape(4,1)
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
        self.selMat = survival / np.sum(survival)
        
        return self.selMat
        



        
        
    def socialMigration(self):
        # matrix here should be selMat
        
        matrix = self.selMat
        
        ms = self.socMigration
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
        
        self.migMat = matrix/np.sum(matrix)

        return self.migMat
        
        
        
        
        
            
                            
    def available(self,matrix,PhiS):
        
        s = self.s
        e = self.e
        
        lm = self.latency[0]
        lf = self.latency[1]
        
        phip = PhiS[0]
        phi = PhiS[1]
        pa1 = np.sum(matrix[:,range(6,8)])/np.sum(matrix[:,range(4,8)]) 
        
        alpha_res = 1/(1+(s*e*phip*lm)/(1-s*lm)) # to replace with mathematica solution for dr
        alpha_a1 = 1/(1 + (s*e*alpha_res*lf)/(1 - s*lf))
        alpha_a0 = 1/(1 + s*e*alpha_res*(1 - phip)*lf/(1 - s*lf))
        alpha_mut = 1/(1 + s*e*lm*(pa1*alpha_a1 + (1 - pa1)*alpha_a0*(1 - phi))/(1-s*lm))
        
        self.availabilities = np.array([alpha_res,alpha_mut,alpha_a0,alpha_a1])
        return self.availabilities
        
        
        
        
                            
        
    def phenoRep(self,matrix,phiVal):
        #here, matrix is the submatrix for one class, one allele of choosiness with associated phi value, i.e. 2 rows of genotypes
        
        availability = self.availabilities
        
        b = self.b
        
        a0 = availability[2]*(1-phiVal)
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
        # matrix here is the matrix for one social class only
        
        phip = self.phiResidentMutant[0]
        phi = self.phiResidentMutant[1]
        
        hi = self.Hattractive * self.Hdominant
        ha = self.Hattractive * (1-self.Hdominant)
        hd = (1-self.Hattractive) * self.Hdominant
        hdelta = (1-self.Hattractive) * (1-self.Hdominant) 
        
        heritp1 = np.array([hi,1/2*(ha+hi),1/2*(hd+hi),1/2*(hdelta+hi),ha,1/2*(hd+ha),1/2*(hdelta+ha),hd,1/2*(hdelta+hd),hdelta])
        heritp2 = np.array([ha,1/2*(hi+ha),1/2*(hdelta+ha),1/2*(hd+ha),hi,1/2*(hdelta+hi),1/2*(hd+hi),hdelta,1/2*(hd+hdelta),hd])
        heritp3 = np.array([hd,1/2*(hdelta+hd),1/2*(hi+hd),1/2*(ha+hd),hdelta,1/2*(hi+hdelta),1/2*(ha+hdelta),hi,1/2*(ha+hi),ha])
        heritp4 = np.array([hdelta,1/2*(hd+hdelta),1/2*(ha+hdelta),1/2*(hi+hdelta),hd,1/2*(ha+hd),1/2*(hi+hd),ha,1/2*(hi+ha),hi])
        
        herit = np.transpose(np.vstack((heritp1,heritp2,heritp3,heritp4)))
        
        phenRep = np.empty([10,4])
        phenRep[:,(0,2)] = self.phenoRep(matrix[(0,2),:],phip)
        phenRep[:,(1,3)] = self.phenoRep(matrix[(1,3),:],phi)
        
        phenNew = np.empty([4,4])
        
        for phen in range(4):
            phenNew[phen,] = np.sum(np.multiply(phenRep,herit[phen,:]),0)
            
        return phenNew

        
        
        
        
        
        
    def matingGenotype(self,submatrix):
        #matrix here is the matrix for 1 social class only
        classFreq = np.sum(submatrix)
        matrix = submatrix/classFreq
        
        phip = self.phiResidentMutant[0]
        phi = self.phiResidentMutant[1]

        s = self.s
        e = self.e
        
        r = self.r

        availability = self.available(matrix,[phip,phi])
        # [res, mut, a0, a1]
        phenGenMat = self.matingPhenotype(matrix)
        
        matFemales = .5*phenGenMat
        matMales = .5*phenGenMat
        
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
        matNew = classFreq * matNew/np.sum(matNew)
        
        return matNew    
        
    
    def linkageDisequilibrium(self):
                
        self.LD = (1-self.r)*self.LD
        
        
        
        
    def downRegulation(self):
        
        repmatrix = self.repMat
        selmatrix = self.selMat
        
        c2Sel = selmatrix[:,range(8,12)]
                          
        MupTot = self.socMigration[0]*(1+self.socMigration[2]/self.socMigration[1])
        HypTot = self.hypergyny[0]*(1+self.hypergyny[2]/self.hypergyny[1])
        
        repmatrix[:,range(4)] = repmatrix[:,range(4)] - c2Sel * (MupTot + HypTot)
        repmatrix[:,range(8,12)] = repmatrix[:,range(8,12)] + c2Sel * (MupTot + HypTot)
               
        self.regMat = repmatrix/np.sum(repmatrix)
        
        return self.regMat
from __future__ import division

NbOfStrategies = 50
mutationRate = 0.001
stepMutation = 1/NbOfStrategies

pUpTown = .1 # proportion of individuals in upper class at t=0
sexRatio = .5 # proportion of males in the population
        
muDiet = [.8,.4] # mu class 1, mu class 2
sigmaDiet = [.2,.4] # sigma class 1, sigma class 2
maxDiet = [1.0,.8] # max survival class 1, max survival class 2

survivalRate = .8
encounterRate = .5
latencyFemale = .8
latencyMale = .6

probaB = 1 # initial frequency of Resident allele for nutrition adaptation
probaC = 1 # initial frequency of Resident allele for choosiness
probad0 = .5 # initial frequency of less socially dominant individuals
probaa0 = .5 # initial frequency of less attractive individuals

# class 1 viability
mu1 = 0.5
sigma1 = 0.1
max1 = 1
# class 2 viability
mu2 = 0.3
sigma2 = 0.2
max2 = 0.8

migrationParamLow = 1
migrationParamHigh = 1.1

pUpTown = .1 # starting frequency of uptown population
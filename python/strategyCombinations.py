# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 14:18:07 2017

@author: Claire
"""

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
homeModules = 'C:\\Users\\Claire\\Documents\\GitHub\\HumanEvol\\python'
sys.path.append(homeModules)
from Parameters import *
import Utilities as ut
from Population import *


dirname = 'C:\\Users\\Claire\\Documents\\GitHub\\HumanEvol'


residentDiets = np.linspace(0,1,NbOfStrategies)
residentPhis = np.linspace(0,1,NbOfStrategies)

allParams = [list(residentDiets),list(residentPhis)]
allParamCombinations = list(it.product(*allParams))

import operator
getX1 = operator.itemgetter(0)
getX2 = operator.itemgetter(1)
x1 = list(map(getX1, allParamCombinations))
x2 = list(map(getX2, allParamCombinations))

with open('%s\\strategies2.txt' % dirname, 'w') as fp:
    for value in range(len(allParamCombinations)):
        fp.writelines('{0},{1}\n'.format(x1[value],x2[value]))
fp.close()



# -*- coding: utf-8 -*-
"""
Created on Tue Aug 01 12:23:34 2017

@author: Claire
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np

def Viability(xVal, muVal, sigmaVal, maxVal = 1):
    
    v = maxVal * np.exp(-(xVal - muVal) ** 2 / (2 * sigmaVal ** 2))
    return v

def ReproductiveSuccess(choosiness):
    f = - choosiness^2 + choosiness + 3/4 
    return(f)
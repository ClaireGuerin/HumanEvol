# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 13:49:04 2017

@author: Claire
"""

import numpy as np

def Quantity(choosiness):
    np.exp(-choosiness)

g1Next = (1-r) * g1 + r * p * q
g2Next = (1-r) * g2 + r * p * (1-q)
g3Next = (1-r) * g3 + r * (1-p) * q
g4Next = (1-r) * g4 + r * (1-p) * (1-q)


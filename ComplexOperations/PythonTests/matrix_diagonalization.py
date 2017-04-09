#!/usr/bin/python

import numpy as np 

#============================

size = 4

A = np.array([[0, 1], [1, 2]])

eigVal, eigVect = np.linalg.eig(A)
print eigVal
print eigVect
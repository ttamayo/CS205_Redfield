#!/usr/bin/python

import matplotlib.pyplot as plt 
import numpy as np 

fileName = 'test'

numbers = []
content = open(fileName, 'rb')
for line in content:
	if not line.startswith('#') and not len(line) < 3:
		numbers.append([float(element) for element in line.split()])
content.close()

array = np.array(numbers)
#print array.shape

for index in range(1, array.shape[1]):
	plt.plot(array[:, 0], array[:, index], lw = 3, label = index - 1)
plt.legend(frameon = False, loc = 'best')
plt.show()
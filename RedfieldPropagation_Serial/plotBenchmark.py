#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns

#==========================================

fileName = 'benchmark.dat'

values = []
content = open(fileName, 'rb')
for line in content:
	if not line.startswith('#'):
		values.append([float(element) for element in line.split()])
content.close()

values = np.array(values)

x = values[:,0][np.where(values[:, 1] != 0)]
y = values[:,1][np.where(values[:, 1] != 0)]

print x, y
print np.polyfit(np.log(x), np.log(y), 1)


for index in range(1, values.shape[1]):
	plt.plot(values[:, 0], values[:, index], marker = 'o')
plt.show()
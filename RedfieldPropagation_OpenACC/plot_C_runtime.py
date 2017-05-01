#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns

sns.set_context('paper', font_scale = 2.0, rc = {'lines.linewidth': 3})

#==========================================

labels = ['C implementation']

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
y1 = values[:,1][np.where(values[:, 1] != 0)]

print 'serial', np.polyfit(np.log(x), np.log(y1), 1)
#print 'parallel', np.polyfit(np.log(x), np.log(y2), 1)


for index in range(1, 2):
	plt.plot(values[:, 0], values[:, index], marker = 'o', label = labels[index - 1])
plt.legend(loc = 'best')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Runtime [s]')
plt.xlabel('Hamiltonian size (n x n)')
plt.savefig('runtimes_C.png', bbox_inches = 'tight')
plt.show()
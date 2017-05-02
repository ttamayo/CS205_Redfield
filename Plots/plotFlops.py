#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns

sns.set_context('paper', font_scale = 2.0, rc = {'lines.linewidth': 3})

#==========================================

labels = ['serial', 'parallel']
colors = ['k', 'r']

#==========================================

fileName = 'benchmark.dat'

values = []
content = open(fileName, 'rb')
for line in content:
	if not line.startswith('#'):
		values.append([float(element.replace(',', '')) for element in line.split()])
content.close()

values = np.array(values)

x = values[:,0][np.where(values[:, 1] != 0)] + 2.
flops = 54 * x**2 + 26 * x**3 + 2 * x**4 + 10*x**6
z_cpu = values[:,1][np.where(values[:, 1] != 0)]
#y_cpu = values[:,3][np.where(values[:, 1] != 0)]
z_gpu = values[:,2][np.where(values[:, 1] != 0)]
#y_gpu = values[:,4][np.where(values[:, 1] != 0)]

plt.plot(x,  flops / z_cpu, label = 'CPU', marker = 'o', markersize = 8)
plt.plot(x,  flops / z_gpu, label = 'GPU', marker = 'o', markersize = 8)

#plt.plot([values[:, 0][0], values[:, 0][-1]], [1.0, 1.0], ls = '--', color = 'k')
#plt.plot(values[:, 0], values[:, 1] / values[:, 2], marker = 'o', markersize = 8, color = 'g')

#for index in range(1, values.shape[1]):
#	plt.plot(values[:, 0], values[:, index], marker = 'o', markersize = 8, label = labels[index - 1], color = colors[index - 1])
plt.legend(loc = 'best')
plt.xlabel('Problem size (n x n)')
plt.ylabel('FLOPs / second')
#plt.xscale('log')
#plt.yscale('log')
#plt.savefig('speedup.png', bbox_inches = 'tight')
plt.show()

#!/usr/bin/python

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
		values.append([float(element) for element in line.split()])
content.close()

values = np.array(values)

x = values[:,0][np.where(values[:, 1] != 0)] + 2.
y1 = values[:,1][np.where(values[:, 1] != 0)]
y2 = values[:,2][np.where(values[:, 1] != 0)]

print 'serial', np.polyfit(np.log(x), np.log(y1), 1)
print 'parallel', np.polyfit(np.log(x), np.log(y2), 1)


plt.plot([values[:, 0][0], values[:, 0][-1]], [1.0, 1.0], ls = '--', color = 'k')
plt.plot(values[:, 0], values[:, 1] / values[:, 2], marker = 'o', markersize = 8, color = '#8A8A8A', label = 'CPU/GPU')

#for index in range(1, values.shape[1]):
#	plt.plot(values[:, 0], values[:, index], marker = 'o', markersize = 8, label = labels[index - 1], color = colors[index - 1])
plt.legend(loc = 'best')
plt.xlabel('Hamiltonian size (n x n)')
plt.ylabel('Speed up')
#plt.xscale('log')
#plt.yscale('log')
plt.savefig('speedup.png', bbox_inches = 'tight')
plt.show()
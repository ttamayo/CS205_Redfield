#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import interpolate

sns.set_context('paper', font_scale = 2.0, rc = {'lines.linewidth': 3})


#values = []
#content = open('test', 'rb')
#for lineIndex, line in enumerate(content):
#	if not line.startswith('#'):
#		values.append([float(element) for element in line.split()])
#content.close()

values = np.loadtxt('test')
print values.shape
#print values[0, 0]
print len(values[0])
print len(values[1])


for i in range(1, values.shape[1]):
	if i == 2 or values.shape[1] - 1:
		ls = '-'
	else:
		ls = '-'
	print i
	plt.plot(values[:, 0], values[:, i],  ls = ls)
plt.plot(values[:, 0], np.sum(values[:, 1:], axis = 1), label = 'Total')

#plt.legend(loc = 'best')
plt.xlabel('Time [fs]')
plt.ylabel('Population')

plt.ylim(-0.05, 1.05)
plt.savefig('RungeKutta_test.png', bbox_inches = 'tight')
plt.show()





#	for i in range(hamiltonian.shape[0]):
#		if 0 < i < hamiltonian.shape[0] - 1:
#			ls = '-'
#		else:
#			ls = '--'
#		plt.plot(timeDomain, np.real(propagator.populations[:, i]), label = i, lw = 3, ls = ls)
#	plt.plot(timeDomain, np.real(np.sum(propagator.populations, axis = 1)), label = 'total')
#	plt.xlabel('Time [fs]')
#	plt.ylabel('Population')
#	plt.legend(ncol = 4, loc = 'best')
#	plt.savefig('population_python.png', bbox_inches = 'tight')
#	plt.show()
#

#	quit()
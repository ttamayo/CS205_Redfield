#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
sns.set_context('paper', font_scale = 2.0)

#======================================================================

name_complex = 'runtimes_mat_mul_complex.dat'
name_complexified = 'runtimes_mat_mul_complexified.dat'
fileNames = [name_complex, name_complexified]
labels = ['complex', 'separated']

#======================================================================

for fileIndex, fileName in enumerate(fileNames):
	content = open(fileName, 'rb')
	for line in content:
		if line.startswith('#'):
			sizes, means, stds = [], [], []
		else:
			linecontent = line.split()
			sizes.append(int(linecontent[0]))
			means.append(np.mean([float(element) for element in linecontent[1:]]))
			stds.append(np.std([float(element) for element in linecontent[1:]]))
	content.close()
	plt.errorbar(sizes, means, yerr = stds, markersize = 10, marker = 'o', lw = 3, label = labels[fileIndex])

plt.xlabel('Matrix size')
plt.ylabel('Runtime [s]')
plt.legend(loc = 'best', frameon = False)
plt.savefig('runtime_comparison_complex_vs_separated.png', bbox_inches = 'tight')
plt.show()

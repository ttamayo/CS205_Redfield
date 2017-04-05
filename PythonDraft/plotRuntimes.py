#!/usr/bin/python

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns 
from colors import colors
sns.set_context('paper', font_scale = 2.0, rc = {'lines.linewidth': 3})

#=====================================================

fileName = 'runtimes.dat'

sizes = []
naive = []
vectorized = []

content = open(fileName, 'rb')
for line in content:
	if not line.startswith('#'):
		linecontent = line.split()
		sizes.append(int(linecontent[0]))
		if linecontent[1] == '-':
			naive.append(np.nan)
		else:
			naive.append(float(linecontent[1]))
		if linecontent[2] == '-':
			vectorized.append(np.nan)
		else:
			vectorized.append(float(linecontent[2]))
content.close()

sizes = np.array(sizes)
naive = np.array(naive)
vectorized = np.array(vectorized)

print np.polyfit(np.log(sizes), np.log(naive), 1)
print np.polyfit(np.log(sizes), np.log(vectorized), 1)

plt.plot(sizes, naive, label = 'naive', marker = 'o', color = colors['darkRed'], markersize = 8)
plt.plot(sizes, vectorized, label = 'vectorized', marker = 'o', color = colors['darkGreen'], markersize = 8)
plt.xlabel('Hamiltonian size (nxn)')
plt.ylabel('Runtime [s]')
plt.xlim(1.9, 16.1)
plt.ylim(0.8 * naive[0], 1.1 * naive[-1])
plt.legend(loc = 'best')
plt.yscale('log')
plt.xscale('log')
plt.savefig('runtimes_loglog.png', bbox_inches = 'tight')
plt.show()

# now we get speed up plots
#plt.clf()

#plt.plot(sizes, naive / vectorized, color = colors['darkBlack'], marker = 'o', markersize = 8)
#plt.xlabel('Hamiltonian size (nxn)')
#plt.ylabel('Speed up')
#plt.savefig('speed_up.png', bbox_inches = 'tight')
#plt.show()
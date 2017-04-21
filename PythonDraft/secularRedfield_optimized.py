#!/usr/bin/python

__author__ = 'Flo'
__email__  = 'fhase@g.harvard.edu'
__date__   = '03/30/17'

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import interpolate

sns.set_context('paper', font_scale = 2.0, rc = {'lines.linewidth': 3})


# for reproducibility
np.random.seed(100691)

#=============================================================================
# unit conversion factors for the code below
# units used in this calculation
# energy       ... cm-1
# time         ... fs
# temperature  ... K

# note: cm1 means cm^{-1}, fs1 means fs^{-1}
s_to_fs = 10.**(15)
eV_to_cm1 = 8065.54429
fs1_to_cm1 = 33356.40952
cm1_to_fs1 = 1 / 33356.40952

#=============================================================================
# physical constants, converted to the above unit system
# TODO: temperature should be a user defined parameter

hbar = 6.582119514 * 10**(-16) 		# eV * s
hbar *= eV_to_cm1 * s_to_fs			# cm-1 * fs
temperature = 300. 					# K
kB = 8.6173303 * 10**(-5) 			# eV / K
kB *= eV_to_cm1
beta = 1 / (kB * temperature)		# 1 / cm-1

#=============================================================================

shape = 8

# initialize hamiltonian as (|0>, |1>, ..., |n>, |RC>)
hamiltonian = np.zeros((shape + 2, shape + 2), dtype = np.complex128)

# we add random excitonic energies and couplings - who cares about the real deal for this example?!
# note: we populate i != j states twice, due to the way the sums are set up
#       not a big deal for now, though
for i in range(1, shape + 1):
	for j in range(1, shape + 1):
		# draw a random number for energy or coupling
		number = np.random.normal(0.0, 500.0)
		hamiltonian[i, j] = number
		hamiltonian[j, i] = number

# now we can define links to targets and losses
# these links are defined in units of 1/fs, as opposed to the hamiltonian above, which is defined in cm1
# target and loss links are defined as dictionaries, with linked excitonic site as index and link rate as target
link_to_target = {-2: 5. / 1000.}  										# we only link to the last excitonic state here, 
																		# but there is no reason to not have more targets
link_to_loss   = { i: 0.5 * 1. / 1000. for i in range(1, shape + 1)}  	# we link all sites to the loss state
																		# but again, this is not universal

#========================================================================================================================

class RedfieldPropagator(object):
	def __init__(self, rho_init, hamiltonian, specParams, timeDomain, link_to_loss, link_to_target):
		self.rho_init       = rho_init		
		self.hamiltonian    = hamiltonian
		self.specParams     = specParams 
		self.timeDomain     = timeDomain
		self.link_to_loss   = link_to_loss
		self.link_to_target = link_to_target

		self.populations       = np.empty((self.timeDomain.shape[0], shape + 2), dtype = np.complex128)
		self.populations[0, :] = np.diag(self.rho_init)

		self._diagonalizeHamiltonian()
		self._getRates()
		self._getRateMatrices()
		self.rho               = np.dot(self.eigVect_inv, np.dot(rho_init, self.eigVect))


	def _diagonalizeHamiltonian(self):
		self.eigVal, self.eigVect = np.linalg.eig(self.hamiltonian)
		self.eigVect_inv = np.linalg.inv(self.eigVect)

	def _getOmegas(self):
		self.omegas = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
		for M in xrange(1, shape + 1):
			for N in xrange(1, shape + 1):
				self.omegas[M, N] = (self.eigVal[M] - self.eigVal[N])

	def _getPhononStatistics(self):
		self.phononStatistics = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
		self.phononStatistics[np.abs(np.real(self.omegas)) > 10**(-12)] = 1 / (np.exp(beta * cm1_to_fs1 * hbar * self.omegas[np.abs(np.real(self.omegas)) > 10**(-12)]) - 1)

	def _getSpectralDensities(self):
		self.spectralDensities = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
		for param_set in self.specParams:
			lambda_k, nu_k, omega_k = param_set[0], 1/param_set[1], param_set[2]
			nu_k *= fs1_to_cm1
			self.spectralDensities += nu_k * lambda_k * self.omegas / (nu_k**2 + (self.omegas - omega_k)**2)
			self.spectralDensities += nu_k * lambda_k * self.omegas / (nu_k**2 + (self.omegas + omega_k)**2)
			self.spectralDensities += np.diag([2 * lambda_k / (nu_k * cm1_to_fs1 * beta * hbar) for i in xrange(shape + 2)])
		self.spectralDensities[0, 0] = 0.
		self.spectralDensities[-1, -1] = 0.

	def _getRates(self):
		self._getOmegas()
		greater = np.where(np.real(self.omegas) > 10**-12)
		smaller = np.where(np.real(self.omegas) < -10**-12)

		self._getPhononStatistics()
		self._getSpectralDensities()
		self.rates = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
		self.rates          = np.copy(self.spectralDensities)
		self.rates[smaller] = (np.transpose(self.spectralDensities) * (np.transpose(self.phononStatistics) + 1.))[smaller]
		self.rates[greater] = (self.spectralDensities * self.phononStatistics)[greater]
		self.rates *= 2 * np.pi * cm1_to_fs1


		for m in self.link_to_loss.keys():
			self.rates[0, m] = self.link_to_loss[m]
		for m in self.link_to_target.keys():
			self.rates[-1, m] = self.link_to_target[m]

		# rates should be ok

	def _get_V(self, i, k):
		matrix = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
		matrix[i, k] += 1.
		matrix_rotated = np.dot(self.eigVect_inv, np.dot(matrix, self.eigVect))
		return matrix_rotated


	def _getRateMatrices(self):
		# rateMatrices are of shape (m, M, N)
		# get them all before going into the propagation
		rm_V        = np.zeros((shape + 2, shape + 2, shape + 2, shape + 2, shape + 2), dtype = np.complex128)
		rm_V_dagg   = np.zeros((shape + 2, shape + 2, shape + 2, shape + 2, shape + 2), dtype = np.complex128)
		rm_V_dagg_V = np.zeros((shape + 2, shape + 2, shape + 2, shape + 2, shape + 2), dtype = np.complex128)

		for m in range(1, shape + 1):
			for M in range(1, shape + 1):
				for N in range(1, shape + 1):
					V = self._get_V(m, m)
					rate = self.rates[M, N]
					rm_V[m, M, N, :, :] = rate * V
					rm_V_dagg[m, M, N, :, :] = np.conj(np.transpose(V))
					rm_V_dagg_V[m, M, N, :, :] = rate * np.dot(np.conj(np.transpose(V)), V)

		for m in self.link_to_target.keys():
			V = self._get_V(-1, m)
			rate = self.link_to_target[m]
			rm_V[-1, m, 0, :, :] = rate * V
			rm_V_dagg[-1, m, 0, :, :] = np.conj(np.transpose(V))
			rm_V_dagg_V[-1, m, 0, :, :] = rate * np.dot(np.conj(np.transpose(V)), V)

		for m in self.link_to_loss.keys():
			V = self._get_V(0, m)
			rate = self.link_to_loss[m]
			rm_V[0, m, 0, :, :] = rate * V
			rm_V_dagg[0, m, 0, :, :] = np.conj(np.transpose(V))
			rm_V_dagg_V[0, m, 0, :, :] = rate * np.dot(np.conj(np.transpose(V)), V)

#		print rm_V
#		quit()

		self.rm_V        = rm_V
		self.rm_V_dagg   = rm_V_dagg
		self.rm_V_dagg_V = rm_V_dagg_V


	def propagate(self):
		dt = self.timeDomain[1] - self.timeDomain[0]
		systemHamiltonian = np.diag(self.eigVal)

		for timeIndex, time in enumerate(self.timeDomain[1:]):

			first = (-1.j / hbar) * (np.dot(systemHamiltonian, self.rho) - np.dot(self.rho, systemHamiltonian))
			
			second = np.zeros((shape + 2, shape + 2), dtype = np.complex128)
			for m in range(1, shape + 1):
				for M in range(1, shape + 1):
					for N in range(1, shape + 1):
						second += np.dot(np.dot(self.rm_V[m, M, N], self.rho), self.rm_V_dagg[m, M, N])
						second -= np.dot(self.rm_V_dagg_V[m, M, N], self.rho) / 2.
						second -= np.dot(self.rho, self.rm_V_dagg_V[m, M, N]) / 2.

			for m in self.link_to_target.keys():
				second += np.dot(np.dot(self.rm_V[-1, m, 0], self.rho), self.rm_V_dagg[-1, m, 0])
				second -= np.dot(self.rm_V_dagg_V[-1, m, 0], self.rho) / 2.
				second -= np.dot(self.rho, self.rm_V_dagg_V[-1, m, 0]) / 2.

			for m in self.link_to_loss.keys():
				second += np.dot(np.dot(self.rm_V[0, m, 0], self.rho), self.rm_V_dagg[0, m, 0])
				second -= np.dot(self.rm_V_dagg_V[0, m, 0], self.rho) / 2.
				second -= np.dot(self.rho, self.rm_V_dagg_V[0, m, 0]) / 2.

			self.rho += (first + second) * dt 

			# recording
			probe = np.dot(np.dot(self.eigVect, self.rho), self.eigVect_inv)
			self.populations[timeIndex + 1, :] = np.diag(probe)

#========================================================================================================================

if __name__ == '__main__':
	timeDomain = np.arange(0., 200., 0.025)
	specParams = [[35., 50., 0.]]
	rho = np.zeros((hamiltonian.shape[0], hamiltonian.shape[0]), dtype = np.complex128)
	rho[1, 1] += 1.

#	print 'hamiltonian'
#	print hamiltonian

	import time
	start = time.time()

	propagator = RedfieldPropagator(rho, hamiltonian, specParams, timeDomain, link_to_loss, link_to_target)
	propagator.propagate()

	print 'Elapsed time: %.5f s' % (time.time() - start)

#	quit()

	for i in range(hamiltonian.shape[0]):
		if 0 < i < hamiltonian.shape[0] - 1:
			ls = '-'
		else:
			ls = '--'
		plt.plot(timeDomain, np.real(propagator.populations[:, i]), label = i, lw = 3, ls = ls)
#	plt.plot(timeDomain, np.real(np.sum(propagator.populations, axis = 1)), label = 'total')
	plt.xlabel('Time [fs]')
	plt.ylabel('Population')
	plt.legend(ncol = 4, loc = 'best')
	plt.savefig('population_python.png', bbox_inches = 'tight')
	plt.show()


	quit()
#!/usr/bin/python

__author__ = 'Flo'
__email__  = 'fhase@g.harvard.edu'
__date__   = '03/30/17'

import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
from scipy import interpolate

# for reproducibility
np.random.seed(100691)
np.set_printoptions(precision = 8)

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
# here we define our system
# we set up the system hamiltonian and define couplings to losses and targets
# note, that we can tune the size of the system hamiltonian
#
# the hamiltonian below will be a matrix of size 'shape + 2'
# 'shape' denotes the number of excitonic sites in the system
# but we also need an additional loss-state and an additional target-states
# 
# the loss state will be denoted as |0>
# excitonic states will be denoted as |1>, |2>, ..., |n>, where 'n' is the number of excitonic sites
# the target state will be denoted as |RC>
#
# the hamiltonian matrix is then set up as (|0>, |1>, ..., |n>, |RC>),
# i.e. the first column is the loss state, the last column is the target state

# number of excitonic sites in the system
shape = 2

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

print hamiltonian
# now we can define links to targets and losses
# these links are defined in units of 1/fs, as opposed to the hamiltonian above, which is defined in cm1
# target and loss links are defined as dictionaries, with linked excitonic site as index and link rate as target
link_to_target = {-2: 5. / 1000.}  										# we only link to the last excitonic state here, 
																		# but there is no reason to not have more targets
link_to_loss   = { i: 0.5 * 1. / 1000. for i in range(1, shape + 1)}  	# we link all sites to the loss state
																		# but again, this is not universal

#========================================================================================================================
# method to diagonalize a matrix
# 		receives a two dimensional, possibly complex valued, numpy array
# 		returns eigenvalues and eigenvectors
# note: eigenvalues are not ordered

def diagonalize(matrix):
	eigVal, eigVect = np.linalg.eig(matrix)
	return eigVal, eigVect #np.array([eigVect[:, i] for i in range(matrix.shape[0])])

#========================================================================================================================
# computes the spectral density as a superposition of drude-lorentz spectral densities defined by parameters in 'params'
#		receives omega, the frequency at which the spectral density is evaluated, and params, explained below
# 		returns the total spectral density at the given frequency omega
# note: params is a list of lists in the form [[lambda_1, 1/nu_1, Omega_1], [lambda_2, 1/nu_2, Omega_2], ..., [lambda_k, 1/nu_k, Omega_k]]
#       every set [lambda, 1/nu, Omega] defines a drude-lorentz spectral density
# 		lambda is given in cm1, 1/nu in fs1 and Omega in cm1, so some unit conversion is necessary in this function
# 			who would have thought carrying units through code was easy?!

def spectral_density(omega, params):
	value = 0
	# loop over all parameter sets
	for param_set in params:
		lambda_k, nu_k, omega_k = param_set[0], 1/param_set[1], param_set[2]
		nu_k *= fs1_to_cm1			# UNIT CONVERSION HERE!!!
		value += 2 * nu_k * lambda_k * omega / (nu_k**2 + (omega - omega_k)**2)
	return value

#========================================================================================================================
# calculate the phonon statistics defined as n(omega) = 1 / (exp(beta * hbar * omega) - 1)
# 		receives omega, the frequency for which to evaluate the phonon statistics, given in cm1
# 		returns the phonon statistics as a dimensionless quantity

def get_phonon_statistics(omega):
	omega *= cm1_to_fs1 	    # UNIT CONVERSION HERE!!!
	return 1 / (np.exp(beta * hbar * omega) - 1)

#========================================================================================================================
# here we calculate the rate based on the spectral density, this is defined in Christoph's 2011 paper, somewhere in the appendix
# 		receives the frequency in cm1 and the spectral density parameter set explained above
# 		returns the coupling rate in fs1
# note: the denominator diverges for small omega, so we need to handle this case explicitly in some limit calculations
#       the chosen cut-off of 10**(-12) is arbitrarily chosen, but Taylor expansion said this cut-off is reasonable

def get_rate(omega, params):
	if omega < -10**(-12):
		value = 2 * np.pi * spectral_density(-omega, params) * (get_phonon_statistics(-omega) + 1)
	elif omega > 10**(-12):
		value = 2 * np.pi * spectral_density(omega, params) * get_phonon_statistics(omega)
	else:
		value = 0
		for param_set in params:
			value += 4 * np.pi * param_set[0] * param_set[1] / (beta * hbar)
	value *= cm1_to_fs1
	return value

#========================================================================================================================
# here we calculate coupling matrices - we couple excitonic states to themselves (get_V), targets (get_V_RC) and losses (get_V_0)
# obviously there is no need to have three separate functions, but it might be more instructive this way
# note: matrix rotations into the excitonic bases are necessary since the propagation is computed in the excitonic basis

def get_V(m, eigVectors):
	matrix = np.zeros((eigVectors.shape[0], eigVectors[0].shape[0]), dtype = np.complex128)
	matrix[m, m] += 1.
	matrix_rotated = np.dot(np.linalg.inv(eigVectors), np.dot(matrix, eigVectors))
	return matrix_rotated

def get_V_0(m, eigVectors):
	matrix = np.zeros((eigVectors.shape[0], eigVectors[0].shape[0]), dtype = np.complex128)
	matrix[0, m] += 1.
	matrix_rotated = np.dot(np.linalg.inv(eigVectors), np.dot(matrix, eigVectors))
	return matrix_rotated

def get_V_RC(m, eigVectors):
	matrix = np.zeros((eigVectors.shape[0], eigVectors[0].shape[0]), dtype = np.complex128)
	matrix[-1, m] += 1.
	matrix_rotated = np.dot(np.linalg.inv(eigVectors), np.dot(matrix, eigVectors))
	return matrix_rotated

#========================================================================================================================
# here we don't propagate the system, but instead calculate the change in the density matrix during the time step dt
# this step is explained (in a quite scattered way) in Christoph's 2011 paper (appendix)
# note: we compute the propagator in the excitonic basis

def propagate(rho, energies, vectors, specParams):
	systemHamiltonian = np.diag(energies)
	# get the Liouville part of the excitonic hamiltonian 
	print 'double checking'
	print np.dot(systemHamiltonian, rho) - np.dot(rho, systemHamiltonian)
	first = np.dot(systemHamiltonian, rho) - np.dot(rho, systemHamiltonian)
	first = first * (-1.j / hbar)

	# this is the contribution of the L_{ex-phon} operator
	# therefore we only loop over the exciton states
	second = np.zeros((vectors.shape[0], vectors.shape[0]), dtype = np.complex128)
	for m in range(1, shape + 1):
		for M in range(1, shape + 1):
			for N in range(1, shape + 1):
				omega = (energies[M] - energies[N])
				V = get_V(m, vectors)
				V_dagg = np.conj(np.transpose(V))
#				print 'V'
#				print V
#				print 'V_dagg'
#				print V_dagg
				new = np.dot(np.dot(V, rho), V_dagg) - np.dot(V_dagg, np.dot(V, rho)) / 2. - np.dot(rho, np.dot(V_dagg, V)) / 2.
				rate = get_rate(omega, specParams)
#				print np.dot(np.dot(V, rho), V_dagg)
#				print - np.dot(V_dagg, np.dot(V, rho)) / 2. - np.dot(rho, np.dot(V_dagg, V)) / 2.
#				print new
#				quit()
				second = second + rate * new



	# now we need to jump into the coupling to RC
	# for this we compute L_{ex-RC}
	third = np.zeros((vectors.shape[0], vectors.shape[0]), dtype = np.complex128)
	for m in link_to_target.keys():
		V = get_V_RC(m, vectors)
		V_dagg = np.conj(np.transpose(V))
		new = np.dot(np.dot(V, rho), V_dagg) - np.dot(V_dagg, np.dot(V, rho)) / 2. - np.dot(rho, np.dot(V_dagg, V)) / 2.
		rate = link_to_target[m]
		third = third + rate * new

	# now we need to work on losses
	fourth = np.zeros((vectors.shape[0], vectors.shape[0]), dtype = np.complex128)
	for m in link_to_loss.keys():
		V = get_V_0(m, vectors)
		V_dagg = np.conj(np.transpose(V))
		new = np.dot(np.dot(V, rho), V_dagg) - np.dot(V_dagg, np.dot(V, rho)) / 2. - np.dot(rho, np.dot(V_dagg, V)) / 2.
		rate = link_to_loss[m]
		fourth = fourth + rate * new

	print 'vectors'
	print vectors
	print 'rho'
	print rho
	print 'first'
	print first
	print 'second'
	print second
	quit()

	return first + second #+ third + fourth

#========================================================================================================================


if __name__ == '__main__':
#	print hamiltonian 			# for reference
	# diagonalize first
	eigEnergies, eigVectors = diagonalize(hamiltonian)
	rho = np.zeros((hamiltonian.shape[0], hamiltonian.shape[0]), dtype = np.complex128)
	# we start population dynamics from site 1 - this is arbitrary
	rho[1, 1] += 1.

	print 'rho\n', rho
	print 'test'
	print np.dot(rho, eigVectors)

	# rotate into excitonic basis
	rho_rotated = np.dot(np.dot(np.linalg.inv(eigVectors), rho), eigVectors)

	# some more parameters
	end = 200.
	dt = 0.025
	timeDomain = np.arange(0, end, dt)
	specDensityParams = [[35., 50., 0.]]
	pop = []
	for i in range(hamiltonian.shape[0]):
		pop.append([rho[i, i]])

	import time
	start = time.time() 

	# here we start the propagation
	for timeIndex in timeDomain[1:]:
		delta_rho = propagate(rho_rotated, eigEnergies, eigVectors, specDensityParams)
		# update according to Euler - this must be changed, Euler is way too inaccurate!!!
		rho_rotated = rho_rotated + delta_rho * dt

		# get density matrix in site basis
		probe = np.dot(np.dot(eigVectors, rho_rotated), np.linalg.inv(eigVectors))
		print '\nprobe:\n'
		print probe
		quit()
		# record populations at this time
		for i in range(hamiltonian.shape[0]):
			pop[i].append(probe[i, i])

	print 'Elapsed time: %.5f s' % (time.time() - start)	

	# plot populations
	pop = np.array(pop)
	for i in range(hamiltonian.shape[0]):
		if 0 < i < hamiltonian.shape[0] - 1:
			ls = '-'
		else:
			ls = '--'
		plt.plot(timeDomain, np.real(pop[i]), label = i, lw = 3, ls = ls)
	plt.plot(timeDomain, np.real(np.sum(pop, axis = 0)), label = 'total')
	plt.xlabel('Time [fs]')
	plt.ylabel('Population')
	plt.legend()
	plt.show()



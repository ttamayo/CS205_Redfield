/* 
 * Calculates transition rates and transition matrices
 * for Redfield propagation of open quantum systems
 * 
 * Contains multiple routines for rate calculations
 * 
 * All functions here are executed before the Redfield propagation
 * ... so only once in the entire computation
 *
 * Methods here assumes Drude-Lorentz spectral densities and parameters for them
 *
 * author: Florian Hase
 *
 */


//********************************************************************//
// loading all modules
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"


//********************************************************************//
// defining physical units
// ... energies     ->  cm^-1
// ... time         ->  fs
// ... temperature  ->  K

// unit conversions
#define cm1_to_fs1  1. / 33356.40952
#define fs1_to_cm1  33356.40952
#define eV_to_cm1   8065.54429
#define s_to_fs     1.e15

// physical constants
#define HBAR        6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV    1. / (6.582119514e-16*eV_to_cm1*s_to_fs)
#define KB          8.6173303e-5*eV_to_cm1
#define PI          3.1415926532897932384626433832
#define T           300.0 // FIXME: temperature should be a parameter
#define BETA        1. / (8.6173303e-5*eV_to_cm1 * 300.0)


//********************************************************************//
// computes spectral density for particular site with particular Drude-Lorentz parameters
// not to be called by the any routine outside of this file
//   j          ... site index
//   omega      ... frequency at which spectral density is to be evaluated
//   params     ... spectral density parameters [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
//   num_params ... number of all parameters

double _spectral_density(int j, double omega, double **params, int num_params) {
	double lambda, nu, omega_shift;
	double spec_dens; 
	// extract spectral density parameters and convert to proper units
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	omega_shift = params[j][2];
	// compute spectral density at given omega
	spec_dens  = nu * lambda * omega / (nu * nu + (omega - omega_shift) * (omega - omega_shift));
	spec_dens += nu * lambda * omega / (nu * nu + (omega + omega_shift) * (omega + omega_shift));
	return spec_dens;
}


//********************************************************************//
// computes the derivative of the spectral density explictly at omega = 0
// expression was derived analytically
// not to be called by the any routine outside of this file
//   j          ... site index
//   omega      ... frequency at which spectral density is to be evaluated
//   params     ... spectral density parameters [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
//   num_params ... number of all parameters

double _spectral_density_derivative(int j, double omega, double **params, int num_params) {
	double lambda, nu, omega_shift;
	// extract spectral density parameters and convert to proper units
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	// compute spectral density derivative
	return 2 * lambda / nu;
}


//********************************************************************//
// computes the phonon statistics for phonons of a particular frequency
// not to be called by the any routine outside of this file
//   omega ... phonon frequency

double _phonon_statistics(double omega) {
	double result;
	// convert units
	omega *= cm1_to_fs1;
	// and get result
	result = 1. / (exp(HBAR * omega * BETA) - 1.);
	return result;
}



//********************************************************************//
// computes the transition rate for one particular site from spectral 
// ... density and phonon statistics
// not to be called by the any routine outside of this file
//   j          ... site index
//   omega      ... frequency at which spectral density is to be evaluated
//   params     ... spectral density parameters [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
//   num_params ... number of all parameters

double _get_rate(int j, double omega, double **params, int num_params) {
	double result;
	// check if omega < 0, omega > 0 or omega = 0
	// in the omega = 0 case we took the limit analytically
	// 1e-12 was found to be close enough to zero for this approximation
	if (omega < -1e-12) {
		// get and return rate for omega < 0
		result = 2 * PI * _spectral_density(j, -omega, params, num_params);
		result *= (_phonon_statistics(-omega) + 1) * cm1_to_fs1;
		return result;
	}
	if (omega > 1e-12) {
		// get and return rate for omega > 0
		result = 2 * PI * _spectral_density(j, omega, params, num_params);
		result *= _phonon_statistics(omega) * cm1_to_fs1;
		return result;
	}
	// get and return rate for omega = 0
	result = 2 * PI * KB * T * HBAR_INV;
	result *= _spectral_density_derivative(j, omega, params, num_params);
	return result;
}



//********************************************************************//
// computes a tensor containing all transition rates for all possible 
// ... transition frequencies in the secular Redfield approximation
//   gammas      ...  three dimensional transition rate tensor
//   params      ...  spectral density parameters [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
//   energies    ...  excitonic eigenenergies (eigenvalues of the hamiltonian)
//   num_params  ...  number of all parameters

void get_rates(double ***gammas, double **params,  double *energies, int SIZE) {
	int unsigned j, M, N;
	double rate;
	// iteration over all sites
	for (j = 0; j < SIZE; j++) {
		// iteration over all eigenstates
		for (M = 0; M < SIZE; M++) {
			// iteration over all eigenstates
			for (N = 0; N < SIZE; N++) {
				// don't assign rate for loss transition
				if (j == 0) {
					gammas[M][N][j] = 0.;
				}
				// don't assign rate for target transition
				if (j == SIZE - 1) {
					gammas[M][N][j] = 0.;
				}
				// assign rate for any intra-exciton transition
				if ((j != 0) && (j != SIZE - 1)) {
					rate = _get_rate(j - 1, energies[M] - energies[N], params, SIZE - 2);
					gammas[M][N][j] = rate;
				}
			}
		}
	}
}


//********************************************************************//
// computes a four dimensional tensor containing all possible transition matrices 
// ... for intra-excitonic transitions and transitions to loss and target states
// equation: \sum\limits_{\omega, M, N} c*_M(m) c_N(m) |M><N| \delta(\hbar\omega - E_M + E_N)
//   V        ... 
//   eigVects ...
//   N        ...

void get_V_matrices(double ****V, double **eigVects, int N) {
	int unsigned i, j, k, l;
	double helper[N][N];
	// make sure all entries are zero (important!)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;
				V[k][j][i][1] = 0.;
				V[k][j][i][2] = 0.;
			}
		}
	}
	// iterate over excitonic states
	for (i = 1; i < N - 1; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;     // encodes loss transitions
				V[k][j][i][1] = 0.;     // encodes intra-excitonic transitions
				V[k][j][i][2] = 0.;     // encodes target transitions
			}
		}
		// compute |m><m| (intra-exciton transition)
		V[i][i][i][1] = 1.;
		// rotate |m><m| into superposition of |M><N| 
		// ... first matrix multiplication
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][1] * eigVects[l][k];
				}
			}
		}
		// ... second matrix multiplication
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][1] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][1] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |0><m| (loss transition)
		V[0][i][i][0] = 1.;
		// rotate |0><m| into superposition of |M><N| 
		// ... first matrix multiplication
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][0] * eigVects[l][k];
				}
			}
		}
		// ... second matrix multiplication 
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][0] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][0] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |RC><m| (target transition)
		V[N - 1][i][i][2] = 1.;
		// rotate |RC><m| into superposition of |M><N| 
		// ... first matrix multiplication
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][2] * eigVects[l][k];
				}
			}
		}
		// ... second matrix multiplication
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][2] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][2] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating
	}
}


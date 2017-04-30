/* 
 * Computes rates for Redfield integration
 * Serial code (only executed once)
 * author: Flo
 *
 */

/**********************************************************************/

// loading modules
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers.h"

/**********************************************************************/

// set unit conversions and physical constants
#define cm1_to_fs1 1. / 33356.40952
#define fs1_to_cm1 33356.40952
#define eV_to_cm1  8065.54429
#define s_to_fs    1.e15

#define HBAR 6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV 1. / (6.582119514e-16*eV_to_cm1*s_to_fs)
#define KB   8.6173303e-5*eV_to_cm1
#define PI   3.1415926532897932384626433832
#define T    300.0 // FIXME: temperature should be a parameter
#define BETA 1. / (8.6173303e-5*eV_to_cm1 * 300.0)


/**********************************************************************/
/***  calculates spectral density for site j  *************************/
double _spectral_density(int j, double omega, double **params, int num_params) {
	// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
	double lambda, nu, omega_shift;
	double spec_dens; 
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	omega_shift = params[j][2];
	spec_dens  = nu * lambda * omega / (nu * nu + (omega - omega_shift) * (omega - omega_shift));
	spec_dens += nu * lambda * omega / (nu * nu + (omega + omega_shift) * (omega + omega_shift));
	return spec_dens;
}

/**********************************************************************/
/***  calculates spectral density derivative at omega=0 for site j  ***/
double _spectral_density_derivative(int j, double omega, double **params, int num_params) {
	// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
	double lambda, nu, omega_shift;
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	return 2 * lambda / nu;
}

/**********************************************************************/
/***  get phonon statistics for particular frequency  *****************/
double _phonon_statistics(double omega) {
	double result;
	omega *= cm1_to_fs1;
	result = 1. / (exp(HBAR * omega * BETA) - 1.);
	return result;
}

/**********************************************************************/
/***  compute rate at given frequency omega for given site j  *********/
double _get_rate(int j, double omega, double **params, int num_params) {
	double result;
	// 
	if (omega < -1e-12) {
		result = 2 * PI * _spectral_density(j, -omega, params, num_params);
		result *= (_phonon_statistics(-omega) + 1) * cm1_to_fs1;
		return result;
	}
	if (omega > 1e-12) {
		result = 2 * PI * _spectral_density(j, omega, params, num_params);
		result *= _phonon_statistics(omega) * cm1_to_fs1;
		return result;
	}
	result = 2 * PI * KB * T * HBAR_INV;
	result *= _spectral_density_derivative(j, omega, params, num_params);
	return result;
}

/**********************************************************************/
/***  compute all transition rates for the given system  **************/

void get_rates(double ***gammas, double **params,  double *energies, int SIZE) {
	// gammas ... matrix in which rates will be stored
	// params ... array with spectral density parameters
	// energies ... eigenenergies 
	// num_params ... number of spectral density parameter sets
	int unsigned j, M, N;
	double rate;
	for (j = 0; j < SIZE; j++) 
		for (M = 0; M < SIZE; M++) 
			for (N = 0; N < SIZE; N++) {
				if (j == 0)          // no extra rate for losses
					gammas[M][N][j] = 0.;
				if (j == SIZE - 1)   // no extra rate for targets
					gammas[M][N][j] = 0.;
				if ((j != 0) && (j != SIZE - 1)) {
					rate = _get_rate(j - 1, energies[M] - energies[N], params, SIZE - 2);
					gammas[M][N][j] = rate;
				}
			}
}

/**********************************************************************/
/***  compute all transition matrices for the given system  ***********/
void get_V_matrices(double ****V, double **eigVects, int N) {
	int unsigned i, j, k, l;
	double helper[N][N];
	// set everything to zero
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;
				V[k][j][i][1] = 0.;
				V[k][j][i][2] = 0.;
			}
		}
	}
	for (i = 1; i < N - 1; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;
				V[k][j][i][1] = 0.;
				V[k][j][i][2] = 0.;
			}
		}
		// compute |m><m|
		V[i][i][i][1] = 1.;
		// rotate
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][1] * eigVects[l][k];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][1] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][1] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |0><m|
		V[0][i][i][0] = 1.;
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][0] * eigVects[l][k];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][0] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][0] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |RC><m|
		V[N - 1][i][i][2] = 1.;
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][2] * eigVects[l][k];
				}
			}
		}
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


/**********************************************************************/
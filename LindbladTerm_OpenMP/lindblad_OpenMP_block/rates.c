
// external modules
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// internal modules 
#include "headers.h"
//#include "matrix_generators.h"

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



// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
double _spectral_density(int j, double omega, double *params, int num_params) {
	double lambda, nu, omega_shift;
	double spec_dens; 
	lambda      = params[3 * j];
	nu          = 1/params[3 * j + 1] * fs1_to_cm1;
	omega_shift = params[3 * j + 2];
	spec_dens  = nu * lambda * omega / (nu * nu + (omega - omega_shift) * (omega - omega_shift));
	spec_dens += nu * lambda * omega / (nu * nu + (omega + omega_shift) * (omega + omega_shift));
//	printf("spec_dens %.5f\n", spec_dens);
	return spec_dens;
}


double _spectral_density_derivative(int j, double omega, double *params, int num_params) {
// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
	double lambda, nu, omega_shift;
	lambda      = params[3 * j];
	nu          = 1/params[3 * j + 1] * fs1_to_cm1;
	return 2 * lambda / nu;
}


double _phonon_statistics(double omega) {
	omega *= cm1_to_fs1;
	double result;
	result = 1. / (exp(HBAR * omega * BETA) - 1.);
	return result;
}



double _get_rate(int j, double omega, double *params, int num_params) {
	double result;
//	printf("omega %.10f\n", omega);
	if (omega < - 1e-12) {
		result = 2 * PI * _spectral_density(j, -omega, params, num_params);
		result *= (_phonon_statistics(-omega) + 1) * cm1_to_fs1;
//		printf("### %.5f %.5f\n", omega, result);
		return result;
	}
	if (omega > 1e-12) {
		result = 2 * PI * _spectral_density(j, omega, params, num_params);
		result *= _phonon_statistics(omega) * cm1_to_fs1;
//		printf("### %.5f %.5f\n", omega, result) ;
		return result;
	}
	// this is else (I hope)
	result = 2 * PI * KB * T * HBAR_INV;
	result *= _spectral_density_derivative(j, omega, params, num_params);
	return result;
}




// gammas ... matrix in which rates will be stored
// params ... array with spectral density parameters
// energies ... eigenenergies 
// num_params ... number of spectral density parameter sets
// Nsites2 ... Nsites + 2
void get_rates(double *gammas, double *params, double *energies, int num_params, int Nsites2) {
	int unsigned j, M, N;
	double rate;

//	for (j = 0; j < Nsites2; j++) {
//		printf("ENERGIES %.10f\n", energies[j]);
//	}

	// get rates for inter-exciton state transitions
	for (j = 0; j < Nsites2; j++) {
		for (M = 0; M < Nsites2; M++) {
			for (N = 0; N < Nsites2; N++) {

				if (j == 0) 
					gammas[M + N * Nsites2 + j * Nsites2*Nsites2] = 0;
				if (j == Nsites2 - 1) 
					gammas[M + N * Nsites2 + j * Nsites2*Nsites2] = 0;

//				if (energies[M] != 0 && energies[N] != 0) {

				rate = _get_rate(j - 1, energies[M] - energies[N], params, num_params);

//				} else {
//					rate = 0;
//				}


//				printf("%d %d %.5f %.5f %.5f \n", M, N, energies[M], energies[N], rate);
				gammas[M + N * Nsites2 + j * Nsites2 * Nsites2] = rate;
			}
		}
	}
}


void get_V_matrices(double *V, double *eigvects, int N) {
//	double* restrict helper;
//	helper = (double *) malloc(sizeof(double) * N*N);
	double helper[N*N];
	double sum;
	int unsigned i, j, k, l;

	#pragma omp parallel for 
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) { 
			for (k = 0; k < N; k++) {
				V[k + j * N + i * N*N + 0 * N*N*N] = 0.;
				V[k + j * N + i * N*N + 1 * N*N*N] = 0.;
				V[k + j * N + i * N*N + 2 * N*N*N] = 0.;
			}
		}
	}


	// V is a tensor of shape SIZE x SIZE x SIZE x 3
	for (i = 1; i < N - 1; i++) {
		// first, we generate |m><m| at position 
		// gen zero matrix
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k + j * N + i * N*N + 0 * N*N*N] = 0.;
				V[k + j * N + i * N*N + 1 * N*N*N] = 0.;
				V[k + j * N + i * N*N + 2 * N*N*N] = 0.;
			}
		}
		// set one position to 1 for |m><m|
		V[i + i * N + i * N*N + 1 * N*N*N] = 1.;
		// now rotate eigvect.T * V * eigvect
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j + k * N] = 0;
				for (l = 0; l < N; l++) {
					helper[j + k * N] += V[j + l * N + i * N*N + N*N*N] * eigvects[l + k * N];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j + k * N + i * N*N + N*N*N] = 0;
				for (l = 0; l < N; l++) {
					V[j + k * N + i * N*N + N*N*N] += eigvects[l + j * N] * helper[l + k * N];
				}
			}
		} // done rotating

	// ---

		// set one position to 1 for |0><m|
		V[0 + i * N + i * N*N + 0 * N*N*N] = 1.;
		// now we rotate again
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j + k * N] = 0;
				for (l = 0; l < N; l++) {
					helper[j + k * N] += V[j + l * N + i * N*N + 0 * N*N*N] * eigvects[l + k * N];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j + k * N + i * N*N] = 0;
				for (l = 0; l < N; l++) {
					V[j + k * N + i * N*N] += eigvects[l + j * N] * helper[l + k * N];
				}
			}
		} // done rotating

	// ---		

		// set one position to 1 for |RC><m|
		V[(N - 1) + i * N + i * N*N + 2 * N*N*N] = 1.;
		// now we rotate again
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j + k * N] = 0;
				for (l = 0; l < N; l++) {
					helper[j + k * N] += V[j + l * N + i * N*N + 2 * N*N*N] * eigvects[l + k * N];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j + k * N + i * N*N + 2 * N*N*N] = 0;
				for (l = 0; l < N; l++) {
					V[j + k * N + i * N*N + 2 * N*N*N] += eigvects[l + j * N] * helper[l + k * N];
				}
			}
		} // done rotating

	} 

//	V[i + j * N + k * N*N + (0, 1, 2) * N*N*N]
}



void get_V(double *V, double *eigvects, int i, int k, int N) {
        gen_zero_matrix_real(V, N);
        V[i + k * N] = 1.;
        // now we need to rotate
        rotate(V, eigvects, N);
}




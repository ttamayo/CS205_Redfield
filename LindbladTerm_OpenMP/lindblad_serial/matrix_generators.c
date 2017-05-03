
#include <stdio.h>
#include <stdlib.h>

#include "headers.h"

/********************************************************************/
/*** Generators for real valued matrices ****************************/
/********************************************************************/

// generates the identity matrix of size N x N
void gen_identity_real(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			if (i == j) 
				A[i + j * N] = 1.;
			else
				A[i + j * N] = 0.;
}


// generates a matrix of size N x N with only zeros
void gen_zero_matrix_real(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			A[j + i * N] = 0.;
}


// generates a matrix of size N x N with only zeros
void gen_one_matrix_real(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			A[i + j * N] = 1.;
}


void gen_random_matrix_real(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            // each component between 0 and 9
            mat[i + j * N] = 10;
}



void gen_random_hamiltonian_real(double *H, int N) {
	gen_zero_matrix_real(H, N);
	int unsigned i, j;
	double number;
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			if (i == j) {
				number = rand()%800 - 400;
			}
			else {
				number = rand()%200 - 100;
			}
			H[i + j * N] = number;
			H[j + i * N] = number;
		}
	}
}

/********************************************************************/
/*** Generators for complex valued matrices *************************/
/********************************************************************/


void gen_identity_complex(double *A_real, double *A_imag, int N) {
	gen_identity_real(A_real, N);
	gen_zero_matrix_real(A_imag, N);
}

void gen_zero_matrix_complex(double *A_real, double *A_imag, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A_real[j + i * N] = 0.;
			A_imag[j + i * N] = 0.;
		}
	}
}


/********************************************************************/
/*** Generators for test cases  *************************************/
/********************************************************************/


void gen_test_hamiltonian(double *A) {
	A[1 + 1 * 4] = -29.79868321;
	A[2 + 2 * 4] = -272.80220217;
	A[1 + 2 * 4] =  452.71157134;
	A[2 + 1 * 4] =  452.71157134;
}


void gen_test_spec_densities(double *params, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		params[3 * i] = 35;
		params[3 * i + 1] = 50;
		params[3 * i + 2] = 0;
	}
}


void gen_test_links(double *links_to_loss, double *links_to_target, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		links_to_loss[i] = 0;
		links_to_target[i] = 0;
	}
	for (i = 1; i < N - 1; i++) {
		links_to_loss[i] = 0.001;
		links_to_target[2] = 0.005;
	}
}

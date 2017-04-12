
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
			A[i + j * N] = 0.;
}


// generates a matrix of size N x N with only zeros
void gen_one_matrix_real(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			A[i + j * N] = 1.;
}



void gen_test_hamiltonian(double *A) {
	A[1 + 1 * 4] = -29.79868321;
	A[2 + 2 * 4] = -272.80220217;
	A[1 + 2 * 4] =  452.71157134;
	A[2 + 1 * 4] =  452.71157134;
}


/********************************************************************/
/*** Generators for complex valued matrices *************************/
/********************************************************************/


void gen_identity_complex(double *A_real, double *A_imag, int N) {
	gen_identity_real(A_real, N);
	gen_zero_matrix_real(A_imag, N);
}

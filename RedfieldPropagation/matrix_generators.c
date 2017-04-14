
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



/********************************************************************/
/*** Generators for complex valued matrices *************************/
/********************************************************************/


void gen_identity_complex(double *A_real, double *A_imag, int N) {
	gen_identity_real(A_real, N);
	gen_zero_matrix_real(A_imag, N);
}

void gen_zero_matrix_complex(double *A_real, double *A_imag, int N) {
	gen_zero_matrix_real(A_real, N);
	gen_zero_matrix_real(A_imag, N);
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


void gen_test_spec_densities(double *params) {
	params[0] = 50;
	params[1] = 35;
	params[2] = 0;
	params[3] = 50;
	params[4] = 35;
	params[5] = 0;
}


void gen_test_links(double *links_to_loss, double *links_to_target, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		links_to_loss[i] = 0.0005;
		links_to_target[1] = 0.005;
	}
}

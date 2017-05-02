/*
 * C program for matrix multiplication 
 * with complex matrices written out as 
 * real part matrices and imaginary part
 * matrices
 * 
 * author: Flo
 *
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "headers.h"
//#include "multiple_matrix_operations.h"

/* ----------------------------------------------------------- */
// operations on matrix and scalars
void matrix_mul_scalar(double *A, double scalar, int N) {
 	int unsigned i, j;
 	for (i = 0; i < N; i++)
 		for (j = 0; j < N; j++)
 			A[j + i * N] *= scalar;
 }


/* ----------------------------------------------------------- */
// matrix operations on real valued matrices

// adding two real valued matrices C = A + B (naive implementation)
void matrix_add_real(double *A, double *B, double *C, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			C[j + i * N] = A[j + i * N] + B[j + i * N];
}

// subtracting two real valued matrices C = A - B (naive implementation)
void matrix_sub_real(double *A, double *B, double *C, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			C[j + i * N] = A[j + i * N] - B[j + i * N];
}

// multiplying two real valued matrices C = AB (naive implementation)
void matrix_mul_real(double *A, double *B, double *C, int N) {
	int unsigned i, j, k;
	double sum;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			sum = 0.;
			for (k = 0; k < N; k++) {
				sum += A[i + k * N] * B[k + j * N];
			}
			C[i + j * N] = sum;
		}
	}
}

/* ----------------------------------------------------------- */


void matrix_add_complex(double *A_real, double *A_imag, double *B_real, double *B_imag, double *C_real, double *C_imag, int N) {
	matrix_add_real(A_real, B_real, C_real, N);
	matrix_add_real(A_imag, B_imag, C_imag, N);
}

/* matrix multiplication for complex matrices C = AB 
 * suppose A = a + ib and B = c + id, then
 * C = (a + ib) (c + id) = (ac - bd) + i (ad + bc)
 * which can be rewritten in the form
 * C = [ (a - b) c + b (c - d)] + i [ b (c - d) + (b + a) d]
 * 
 * disadvantage: we need to generate matrices a - b, a + b and c - d
 * 
 * advantage: we only need to perform three instead of four matrix multiplications
 *
 * --> we expect a speed up for large matrices
 *
 */
//#pragma acc routine gang
void matrix_mul_complexified(double *A_real, double *A_imag, double *B_real, double *B_imag, double *C_real, double *C_imag, int N) {

	double *helper_matrix;
	helper_matrix = (double *) malloc(sizeof(double) * (N * N));
	double *h1, *h2, *h3;
	h1 = (double *) malloc(sizeof(double) * (N * N));
	h2 = (double *) malloc(sizeof(double) * (N * N));
	h3 = (double *) malloc(sizeof(double) * (N * N));



	// calculate a - b
	matrix_sub_real(A_real, A_imag, h1, N);
	// calculate c - d
	matrix_sub_real(B_real, B_imag, helper_matrix, N);
	// calculate a + b
	matrix_add_real(A_real, A_imag, h2, N);

	// now we do the multiplications
	// calculate (a - b) c
	matrix_mul_real(h1, B_real, C_real, N);
	// calculate b (c - d)
	matrix_mul_real(A_imag, helper_matrix, h3, N);
	// calculate (b + a) d
	matrix_mul_real(h2, B_imag, C_imag, N);

	// now we calculate real and imaginary part of C
	matrix_add_real(C_real, h3, C_real, N);
	matrix_add_real(C_imag, h3, C_imag, N);

	free((void*) helper_matrix);
	free((void*) h1);
	free((void*) h2);
	free((void*) h3);

}

/* ----------------------------------------------------------- */



	// calculate a - b
//	matrix_sub_real(A_real, A_imag, C_real, N);
	// calculate c - d
//	matrix_sub_real(B_real, B_imag, helper_matrix, N);
	// calculate a + b
//	matrix_add_real(A_real, A_imag, C_imag, N);

	// now we do the multiplications
	// calculate (a - b) c
//	matrix_mul_real(C_real, B_real, C_real, N);
	// calculate b (c - d)
//	matrix_mul_real(A_imag, helper_matrix, helper_matrix, N);
	// calculate (b + a) d
//	matrix_mul_real(C_imag, B_imag, C_imag, N);

	// now we calculate real and imaginary part of C
//	matrix_add_real(C_real, helper_matrix, C_real, N);
//	matrix_add_real(C_imag, helper_matrix, C_imag, N);

//	free((void*) helper_matrix);

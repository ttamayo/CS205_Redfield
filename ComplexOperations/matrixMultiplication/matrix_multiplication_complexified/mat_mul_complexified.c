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

// matrices are of size SIZE*SIZE
#define SIZE 1024

/* ----------------------------------------------------------- */
// methods for initializing matrices 

// generates the identity matrix of size N x N
void get_identity_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			if (i == j)
				A[i + j * N] = 1.;
			else
				A[i + j * N] = 0.;
}

// generate a matrix of size N x N of only zeros
void get_zero_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			A[i + j * N] = 0.;
}

// generate a matrix of size N x N of only ones
void get_one_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			A[i + j * N] = 1.;
}

// print matrix
void print_matrix(double *A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.1f ", A[i + j * N]);
        }
        printf("\n");
    }
    printf("\n");
}


/* ----------------------------------------------------------- */
// matrix operations on real valued matrices

// adding two real valued matrices C = A + B (naive implementation)
void matrix_add_real(double *A, double *B, double *C, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			C[i + j * N] = A[i + j * N] + B[i + j * N];
}

// subtracting two real valued matrices C = A - B (naive implementation)
void matrix_sub_real(double *A, double *B, double *C, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			C[i + j * N] = A[i + j * N] - B[i + j * N];
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
void matrix_mul_complexified(double *A_real, double *A_imag, double *B_real, double *B_imag, double *C_real, double *C_imag, int N) {
	double *helper_matrix;
	helper_matrix = (double *) malloc(sizeof(double) * (N * N));

	// calculate a - b
	matrix_sub_real(A_real, A_imag, C_real, N);
	// calculate c - d
	matrix_sub_real(B_real, B_imag, helper_matrix, N);
	// calculate a + b
	matrix_add_real(A_real, A_imag, C_imag, N);

	// now we do the multiplications
	// calculate (a - b) c
	matrix_mul_real(C_real, B_real, C_real, N);
	// calculate b (c - d)
	matrix_mul_real(A_imag, helper_matrix, helper_matrix, N);
	// calculate (b + a) d
	matrix_mul_real(C_imag, B_imag, C_imag, N);

	// now we calculate real and imaginary part of C
	matrix_add_real(C_real, helper_matrix, C_real, N);
	matrix_add_real(C_imag, helper_matrix, C_imag, N);

	free((void*) helper_matrix);
}

/* ----------------------------------------------------------- */


int main() {
	double *A_real, *A_imag;
	double *B_real, *B_imag;
	double *C_real, *C_imag;
	A_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	A_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	B_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	B_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	C_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	C_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));

	// initialize matrices
	get_identity_matrix(A_real, SIZE);
	get_zero_matrix(A_imag, SIZE);
	get_one_matrix(B_real, SIZE);
	get_zero_matrix(B_imag, SIZE);
	get_identity_matrix(C_real, SIZE);
	get_zero_matrix(C_imag, SIZE);


	clock_t tic = clock();

	// multiply matrices
	matrix_mul_complexified(A_real, A_imag, B_real, B_imag, C_real, C_imag, SIZE);

	clock_t toc = clock();

	double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
	printf("\nTime Elapsed: %f seconds\n", time_spent);
    printf("\n*********************************************\n\n");

    return 0;
}
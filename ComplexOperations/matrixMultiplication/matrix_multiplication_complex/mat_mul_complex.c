/*
 * C program for matrix multiplication 
 * using the standard complex.h library
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
void get_identity_matrix(double complex *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			if (i == j)
				A[i + j * N] = 1.;
			else
				A[i + j * N] = 0.;
}

// generates a matrix of size N x N with all entries set to 1. + 0.i
void get_one_matrix(double complex *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			A[i + j * N] = 1.;
}

// prints a complex matrix
void print_matrix(double complex *A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.1f + i%.1f ", creal(A[i + j * N]), cimag(A[i + j * N]));
        }
        printf("\n");
    }
    printf("\n");
}


/* ----------------------------------------------------------- */

// naive complex matrix multiplication
void matrix_mul(double complex *A, double complex *B, double complex *C, int N) {
	int unsigned i, j, k;
	double complex sum; 

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			sum = 0;
			for (k = 0; k < N; k++) {
				sum += A[i + k * N] * B[k + j * N];
			}
			C[i + j * N] = sum;
		}
	}
}

/* ----------------------------------------------------------- */


int main() {
	double complex *A, *B, *C;
	A = (complex double *) malloc(sizeof(complex double) * (SIZE * SIZE));
	B = (complex double *) malloc(sizeof(complex double) * (SIZE * SIZE));
	C = (complex double *) malloc(sizeof(complex double) * (SIZE * SIZE));

	// initialize matrices
	get_identity_matrix(A, SIZE);
	get_one_matrix(B, SIZE);
	get_identity_matrix(C, SIZE);

	clock_t tic = clock();

	// multiply matrices
	matrix_mul(A, B, C, SIZE);

	clock_t toc = clock();
	double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
	printf("\nTime Elapsed: %f seconds\n", time_spent);
    printf("\n*********************************************\n\n");

    return 0;
}
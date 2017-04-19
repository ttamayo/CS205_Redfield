#include <stdio.h>
#include <stdlib.h>

#include "headers.h"

<<<<<<< HEAD
=======
<<<<<<< HEAD
=======

>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
// load diagonalization routine from lapack
// this method diagonalizes double precision real valued symmetric matrices
// reference: http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/dsyev.html
extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, 
				  double* w, double* work, int* lwork, int* info);


/********************************************************************/


<<<<<<< HEAD
//#pragma acc routine worker
=======
<<<<<<< HEAD
#pragma acc routine worker
void transpose(double *A, int N) {
	int unsigned i, j, index;
	double element;

	#pragma acc data present_or_copy(A[0:N*N]) copyin(element)
	#pragma acc loop independent
	for (i = 0; i < N; i++) {
		#pragma acc loop independent
		for (j = 0; j < N; j++) {
			if (i <= j) {
				element = A[i + j * N];
				A[i + j * N] = A[j + i * N];
				A[j + i * N] = element;
			}
=======
#pragma acc routine
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
void transpose(double *A, int N) {
	int unsigned i, j, index;
	double element;

//	#pragma acc data present_or_copy(A[0:N*N]) copyin(element)
//	#pragma acc loop independent
	for (i = 0; i < N; i++) {
<<<<<<< HEAD
//		#pragma acc loop independent
		for (j = 0; j < N; j++) {
			if (i <= j) {
				element = A[i + j * N];
				A[i + j * N] = A[j + i * N];
				A[j + i * N] = element;
			}
=======
		for (j = i; j < N; j++) {
			element = A[i + j * N];
			A[i + j * N] = A[j + i * N];
			A[j + i * N] = element;
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
		}
	}
}


<<<<<<< HEAD
//#pragma acc routine worker
void rotate(double *A, double *eigvect, int N) {
	double *h;
	h = (double *) malloc(sizeof(double) * N * N);

//	#pragma acc data present_or_copy(A[0:N*N], eigvect[0:N*N]) create(h[0:N*N])
=======
<<<<<<< HEAD
#pragma acc routine worker
void rotate(double *A, double *eigvect, int N) {
	double *h;
	h = (double *) malloc(sizeof(double) * N * N);

	#pragma acc data present_or_copy(A[0:N*N], eigvect[0:N*N]) create(h[0:N*N])
=======
#pragma acc routine
void rotate(double *A, double *eigvect, int N) {
	double *h;
	h = (double *) malloc(sizeof(double) * N * N);
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	matrix_mul_real(A, eigvect, h, N);
	transpose(eigvect, N);
	matrix_mul_real(eigvect, h, A, N);
	transpose(eigvect, N);
<<<<<<< HEAD
//	#pragma end data
=======
<<<<<<< HEAD
	#pragma end data
=======
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

	free((void*) h);
}



void diagonalize(double *A, double *D, int N) {
	int n = N, info, lwork;
	double wkopt;
	double* work;
	int LDA = N;

	// we need to get the optimal work load first
	lwork = -1;
	dsyev_("Vectors", "Upper", &n, A, &LDA, D, &wkopt, &lwork, &info);
	// then we perform the actual diagonalization
	lwork = (int)wkopt;
	work = (double*) malloc(lwork * sizeof(double));
	dsyev_("Vectors", "Upper", &n, A, &LDA, D, work, &lwork, &info);
//	free((void*) work);
}

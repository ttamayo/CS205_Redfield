/*
 * Utility functions for computing the propagation of an open quantum system
 * under the secular Redfield equations
 *
 * none of these functions are performance critical
 * most of these functions contain numbers which should be provided by the user
 *
 * author: Florian Hase
 * 
 */

//********************************************************************//
// loading all modules
#include <stdio.h>
#include <stdlib.h>
#include "headers.h"

//********************************************************************//
// load lapack's diagonalization routine

extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                                  double* w, double* work, int* lwork, int* info);

//********************************************************************//
// diagonalize a given matrix
// note: we implemented matrices as two dimensional objects
//       but lapack requires one dimension arrays, so we need to convert first

void diagonalize(double **matrix, double *diagonal, int N) {
	int n = N, info, lwork;
	double wkopt;
	double* work;
	int LDA = N;

	int unsigned i, j;
	double *flat_matrix;
	flat_matrix = (double *) malloc(sizeof(double) * N*N);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			flat_matrix[i + j * N] = matrix[i][j];
		}
	}

    // we need to get the optimal work load first
    lwork = -1;
    dsyev_("Vectors", "Upper", &n, flat_matrix, &LDA, diagonal, &wkopt, &lwork, &info);
    // then we perform the actual diagonalization
    lwork = (int)wkopt;
    work = (double*) malloc(lwork * sizeof(double));
    dsyev_("Vectors", "Upper", &n, flat_matrix, &LDA, diagonal, work, &lwork, &info);

    for (i = 0; i < N; i++) {
    	for (j = 0; j < N; j++) {
    		matrix[i][j] = flat_matrix[i + j * N];
    	}
    }
}


//********************************************************************//
// generate a matrix with only zeros as entries 

void gen_zero_matrix(double **A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A[i][j] = 0.;
		}
	}
}


//********************************************************************//
// generate a random hamiltonian

void gen_random_hamiltonian_real(double **H, int N) {
	int unsigned i, j;
	double number;
	gen_zero_matrix(H, N);
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			if (i == j) 
				number = rand()%800 - 400; 		// excited state energies
			else
				number = rand()%200 - 100; 		// couplings
			H[i][j] = number;
			H[j][i] = number;
		}
	}
}


//********************************************************************//
// generate loss and target transition rates


void gen_test_links(double *links_to_loss, double *links_to_target, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		links_to_loss[i] = 0.;
		links_to_target[i] = 0.;
	}
	for (i = 1; i < N - 1; i++) {
		links_to_loss[i] = 0.001;
		links_to_target[3] = 0.005;
	}
}


//********************************************************************//
// generate spectral density parameters

void gen_test_spec_densities(double **params, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		params[i][0] = 35.;
		params[i][1] = 50.;
		params[i][2] = 0.;
	}
}

//********************************************************************//
// transpose a matrix

void transpose(double **matrix, int N) {
	int unsigned i, j;
	double element;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i <= j) {
				element = matrix[i][j];
				matrix[i][j] = matrix[j][i];
				matrix[j][i] = element;
			}
		}
	}
}


//********************************************************************//
// rotate a matrix with orthogonal basis

void rotate(double **matrix, double **eigVects, int N) {
	int unsigned i, j, k;
	double helper[N][N];
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			helper[i][j] = 0;
			for (k = 0; k < N; k++) {
				helper[i][j] += matrix[i][k] * eigVects[k][j];
			}
		}
	}
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			matrix[i][j] = 0.;
			for (k = 0; k < N; k++) {
				matrix[i][j] += eigVects[k][i] * helper[k][j];
			}
		}
	}
}



//********************************************************************//
// print a matrix

void print_matrix_real(double **A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.5f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

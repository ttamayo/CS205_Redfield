
#include <stdio.h>
#include <stdlib.h>


// load diagonalization routine from lapack
// this method diagonalizes double precision real valued symmetric matrices
// reference: http://physics.oregonstate.edu/~landaur/nacphy/lapack/routines/dsyev.html
extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, 
				  double* w, double* work, int* lwork, int* info);

/* --------------------------------------------------------------------- */

#define SIZE 2

/* --------------------------------------------------------------------- */

// generates a matrix of size N x N with all entries set to 1.
void get_one_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			A[i + j * N] = 1.;
}

// method to print matrices
void print_matrix(double *A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.3f  ", A[i + j * N]);
        }
        printf("\n");
    }
    printf("\n");
}


/* --------------------------------------------------------------------- */

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
	dsyev_(&JOBZ, &UPLO, &n, A, &LDA, D, work, &lwork, &info);
	free((void*) work);
}

/* --------------------------------------------------------------------- */

// run diagonalization
int main(){
	double *A;
	A = (double *) malloc(sizeof(double) * (SIZE * SIZE));
	double *D;
	D = (double *) malloc(sizeof(double) * SIZE);

	get_one_matrix(A, SIZE);
	int unsigned i;
	for (i = 0; i < SIZE; i++) {
		D[i] = 1.0;
	}
	print_matrix(A, SIZE);
	printf("\n");
	for (i = 0; i < SIZE; i++) {
		printf("%.3f ", D[i]);
	}
	printf("\n");
	printf("-- setup complete --\n");

	diagonalize(A, D, SIZE);
	print_matrix(A, SIZE);
	for (i = 0; i < SIZE; i++) {
		printf("%.3f ", D[i]);
	}
	printf("\n");
	return 0;
}


#include <stdio.h>
#include <time.h>

#define SIZE 512


void get_identity_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			if (i == j)
				A[i + j * N] = 1.;
			else
				A[i + j * N] = 0.;
}


void get_one_matrix(double *A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) 
		for (j = 0; j < N; j++) 
			A[i + j * N] = 1.;
}


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


/* ============================================================= */


void matrix_mul(double *A, double *B, double *C, int N) {
	int unsigned i, j, k;
	double sum; 

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


double cumulative_difference(double *A, double *B, int N) {
	int unsigned i, j;
	double result = 0;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			result += A[i + j * N] - B[i + j * N];
		}
	}
	return result;
}








int main() {
	double A[SIZE * SIZE];
	double B[SIZE * SIZE];
	double C[SIZE * SIZE];

	// initialize matrices
	get_identity_matrix(A, SIZE);
	get_one_matrix(B, SIZE);
	get_identity_matrix(C, SIZE);

//	print_matrix(A, SIZE);
//	print_matrix(B, SIZE);
//	print_matrix(C, SIZE);

	clock_t tic = clock();

	// multiply matrices
	matrix_mul(A, B, C, SIZE);

	clock_t toc = clock();
	double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
	double diff = 0;
    diff = cumulative_difference(B, C, SIZE);
    printf("\ndifference: %f", diff);
	printf("\nTime Elapsed: %f seconds\n", time_spent);
    printf("\n*********************************************\n\n");

    return 0;
}
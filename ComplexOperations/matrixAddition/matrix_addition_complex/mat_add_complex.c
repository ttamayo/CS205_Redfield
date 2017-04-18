/*
 * C program for matrix addition: A + B = C
 * using stanard library for complex numbers
 *
 * author: hannahsim
 *
 */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// matrix size: SIZE * SIZE
#define SIZE 4


// For testing purposes
void get_one_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            mat[i + j * N] = 1. + 1.*I;
}

void get_zero_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0. + 0.*I;
}

double complex sumMatrix(double complex *mat, int N) {
    double complex tot = 0;
    for (int i = 0; i < N * N; ++i) {
        tot += mat[i];
    }
    return tot;
}





// Matrix addition
void matrix_add_complex(double complex *A, double complex *B, double complex *C, int N) {

    int unsigned i,j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++)
            C[i+j*N] = A[i+j*N] + B[i+j*N];
}






// PROGRAM MAIN
int main() {

    // Allocate matrices
    double complex *A;
    double complex *B;
    double complex *C;
    A = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    B = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    C = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));

    // Initialize matrices
    get_one_matrix(A, SIZE);
    get_one_matrix(B, SIZE);
    get_zero_matrix(C, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix addition
    matrix_add_complex(A, B, C, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);


    // Compute expected results and check
    double complex sumC = sumMatrix(C, SIZE);
  
    printf("CHECK:\n");
    printf("------\n");
    printf("Expected sum of either: %f\n", 2.*SIZE*SIZE);
    printf("Sum of real components of C: %f\n",creal(sumC));
    printf("Sum of imag components of C: %f\n",cimag(sumC));
    printf("\n\n");


    return 0;

}




















/*
 * C program for matrix addition: A + B = C
 * Treating real and complex components as separate matrices
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


// Helper Function: Addition operation on real valued matrices
// (original author: flo)
void matrix_add_real(double *A, double *B, double *C, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            C[i + j * N] = A[i + j * N] + B[i + j * N];
}




// For testing purposes
void get_one_matrix(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++) 
            mat[i + j * N] = 1.;
}

void get_zero_matrix(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0.;
}

double sumMatrix(double *mat, int N) {
    double realSum = 0;
    for (int i = 0; i < N * N; ++i) {
        realSum += mat[i];
    }
    return realSum;
}





// Matrix addition: treating real and imaginary components separately
void matrix_add_complexified(double *A_real, double *A_imag, double *B_real, double *B_imag, double *C_real, double *C_imag, int N) {

    // Add real and imaginary parts 
    matrix_add_real(A_real, B_real, C_real, N);
    matrix_add_real(A_imag, B_imag, C_imag, N);
}




// PROGRAM MAIN
int main() {

    // Allocate matrices
    double *AA_real, *AA_imag;
    double *BB_real, *BB_imag;
    double *CC_real, *CC_imag;
    AA_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    AA_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    BB_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    BB_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    CC_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    CC_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));

    // Initialize matrices
    get_one_matrix(AA_real, SIZE);
    get_one_matrix(AA_imag, SIZE);
    get_one_matrix(BB_real, SIZE);
    get_one_matrix(BB_imag, SIZE);
    get_zero_matrix(CC_real, SIZE);
    get_zero_matrix(CC_imag, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix addition
    matrix_add_complexified(AA_real, AA_imag, BB_real, BB_imag, CC_real, CC_imag, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);


    // Compute expected results and check
    double sumCReal = sumMatrix(CC_real, SIZE);
    double sumCImag = sumMatrix(CC_imag, SIZE);
  
    printf("CHECK:\n");
    printf("------\n");
    printf("Expected sum of either: %f\n", 2.*SIZE*SIZE);
    printf("Sum of real components of C: %f\n",sumCReal);
    printf("Sum of imag components of C: %f\n",sumCImag);
    printf("\n\n");


    return 0;

}



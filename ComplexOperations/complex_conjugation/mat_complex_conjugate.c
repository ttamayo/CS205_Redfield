/*
 * C program for matrix complex conjugation:
 * element x+yi --> x-yi
 * using standard library for complex numbers
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





// Matrix complex conjugation
void matrix_complex_conjugation(double complex *A, double complex *Astar, int N) {

    int unsigned i,j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++)
            Astar[i+j*N] = conj(A[i+j*N]);
}






// PROGRAM MAIN
int main() {

    // Allocate matrices
    double complex *A;
    double complex *Astar;
    A = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    Astar = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));

    // Initialize matrices
    get_one_matrix(A, SIZE);
    get_zero_matrix(Astar, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix addition
    matrix_complex_conjugation(A, Astar, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);


    // Compute expected results and check
    double complex sumAstar = sumMatrix(Astar, SIZE);
  
    printf("CHECK:\n");
    printf("------\n");
    printf("Expected sum of real: %i\n", SIZE*SIZE);
    printf("Expected sum of imag: %i\n", -SIZE*SIZE);
    printf("Sum of real components of C: %f\n",creal(sumAstar));
    printf("Sum of imag components of C: %f\n",cimag(sumAstar));
    printf("\n\n");


    return 0;

}

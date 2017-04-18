/*
 * C program for matrix conjugate tranpose
 * element x+yi --> y-xi
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
void get_random_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = rand()%10 + rand()%10*I; // each component between 0 and 9
}

void get_zero_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0. + 0.*I;
}


// Matrix conjugate tranpose
void matrix_conjugate_tranpose(double complex *A, double complex *Adag, int N) {
    int unsigned i,j;
    for (i=0; i<N; ++i)
        for (j=0; j<N; ++j)
            Adag[j+i*N] = conj(A[i+j*N]);
}


// PROGRAM MAIN
int main() {
    
    // Allocate matrices
    double complex *A;
    double complex *Adag;
    A = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    Adag = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));

    // Initialize matrices
    get_random_matrix(A, SIZE);
    get_zero_matrix(Adag, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix conjugate tranposing
    matrix_conjugate_tranpose(A, Adag, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);


/*
    // Check results for small sized matrices
    printf("\n\nMatrix A\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(A[i]),cimag(A[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    printf("\n\nMatrix A^{dag}\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(Adag[i]),cimag(Adag[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }
*/

    
    return 0;

}

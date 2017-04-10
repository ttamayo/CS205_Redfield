/*
 * C program for matrix transpose: 
 * A^T_{ij} = A_{ji}
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





// Matrix tranpose
void matrix_transpose_complex(double complex *A, double complex *AT, int N) {
    int unsigned i,j;
    for (i = 0; i < N; i++) 
        for (j = 0; j < N; j++)
            AT[j+i*N] = A[i+j*N];
}




// PROGRAM MAIN
int main() {

    // Allocate matrices
    double complex *A;
    double complex *AT;
    A = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    AT = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));

    // Initialize matrices
    get_random_matrix(A, SIZE);
    get_zero_matrix(AT, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix addition
    matrix_transpose_complex(A, AT, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);

/*
    // Check results for small matrices
    printf("\n\nMatrix A\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(A[i]),cimag(A[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    } 

    printf("\n\nMatrix A^T\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(AT[i]),cimag(AT[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }
*/

    return 0;

}




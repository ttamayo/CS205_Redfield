/*
 * C program for matrix conjugate tranpose
 * element x+yi --> y-xi
 * by treating real and imag parts separately
 *
 * author: hannahsim
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// matrix size: SIZE * SIZE
#define SIZE 4


// For testing purposes
void get_random_matrix(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = rand()%10; // between 0 and 9
}

void get_zero_matrix(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0.;
}


void matrix_conjugate_transpose_real(double *A_real, double *Adag_real, int N) {
    int unsigned i,j;
    for (i=0; i<N; ++i)
        for (j=0; j<N; ++j)
            Adag_real[j+i*N] = A_real[i+j*N];
}

void matrix_conjugate_transpose_imag(double *A_imag, double *Adag_imag, int N) {
    int unsigned i,j;
    for (i=0; i<N; ++i)
        for (j=0; j<N; ++j)
            Adag_imag[j+i*N] = -A_imag[i+j*N];
}


// Matrix conjugate tranpose
void matrix_conjugate_tranpose_complexified(double *A_real, double *A_imag, double *Adag_real, double *Adag_imag, int N) {
    
    matrix_conjugate_transpose_real(A_real, Adag_real, N);
    matrix_conjugate_transpose_imag(A_imag, Adag_imag, N);
}


// PROGRAM MAIN
int main() {
    
    // Allocate matrices
    double *A_real; 
    double *A_imag;
    double *Adag_real;
    double *Adag_imag;
    A_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    A_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    Adag_real = (double *) malloc(sizeof(double) * (SIZE * SIZE));
    Adag_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));

    // Initialize matrices
    get_random_matrix(A_real, SIZE);
    get_random_matrix(A_imag, SIZE);
    get_zero_matrix(Adag_real, SIZE);
    get_zero_matrix(Adag_imag, SIZE);

    // For timing
    clock_t tic, toc;

    tic = clock();

    // Perform matrix conjugate tranposing
    matrix_conjugate_tranpose_complexified(A_real, A_imag, Adag_real, Adag_imag, SIZE);

    toc = clock();

    // Compute and display elapsed time
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);



    // Check results for small sized matrices
    printf("\n\nMatrix A\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", A_real[i],A_imag[i]);
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    printf("\n\nMatrix A^{dag}\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", Adag_real[i],Adag_imag[i]);
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    
    return 0;

}

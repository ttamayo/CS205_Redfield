/*
 * C program to convert "complexified" to "complex" and vice versa.
 * 
 * Types:
 * "Complexified": using two real matrices for real and imag parts 
 * "Complex"     : using standard complex type matrix
 *
 */ 

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// matrix size: SIZE * SIZE
#define SIZE 4

void get_random_complex_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = rand()%10 + rand()%10*I; // each component between 0 and 9
}

void get_zero_real_matrix(double *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0.;
}


// COMPLEX to COMPLEXIFIED conversion
void complex2complexified(double complex *mat, double *mat_real, double *mat_imag, int N ) {
    int unsigned i,j;
    for (i=0;i<N;++i) {
        for (j=0;j<N;++j) {
            mat_real[i+j*N] = creal(mat[i+j*N]);
            mat_imag[i+j*N] = cimag(mat[i+j*N]);
        }
    }
}

// COMPLEXIFIED to COMPLEX conversion
void complexified2complex(double *mat_real, double *mat_imag, double complex *mat, int N) {
    int unsigned i,j;
    for (i=0;i<N;++i)
        for (j=0;j<N;++j)
            mat[i+j*N] = mat_real[i+j*N] + mat_imag[i+j*N]*I;
}


// PROGRAM MAIN
int main() {

    // Allocate matrices
    double complex *mat;
    double *mat_real;
    double *mat_imag;
    mat = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    mat_real = (double *) malloc(sizeof(double) * (SIZE * SIZE)); 
    mat_imag = (double *) malloc(sizeof(double) * (SIZE * SIZE));


    // Initialize matrices
    get_random_complex_matrix(mat, SIZE);
    get_zero_real_matrix(mat_real, SIZE);
    get_zero_real_matrix(mat_imag, SIZE);

    printf("\n\nINITIAL: Complex Matrix:\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(mat[i]),cimag(mat[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }


    complex2complexified(mat,mat_real,mat_imag,SIZE);
    

    printf("\n\nReal part of Matrix:\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f ",mat_real[i]);
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    printf("\n\nImag part of Matrix:\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f ",mat_imag[i]);
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    get_random_complex_matrix(mat,SIZE); // "shuffle"

    complexified2complex(mat_real,mat_imag,mat,SIZE);

    printf("\n\nFINAL: Complex Matrix:\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%.f + %.fi ", creal(mat[i]),cimag(mat[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }

    return 0;
    
}

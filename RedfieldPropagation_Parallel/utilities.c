
// external modules 
#include <stdio.h>
#include <stdlib.h>

// internal modules
#include "headers.h"


// print matrix
//#pragma acc routine
void print_matrix_real(double *A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.3f ", A[i + j * N]);
        }
        printf("\n");
    }
    printf("\n");
//    printf("done printing\n");
}

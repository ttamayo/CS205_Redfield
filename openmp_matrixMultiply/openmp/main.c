// Matrix-matrix multiplication using OpenMP
// Block implementation

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define SIZE 4096
#define BLOCK_SIZE 16

int main(void) {


	double A[SIZE][SIZE];
	double B[SIZE][SIZE];
	double C[SIZE][SIZE];
	
	int unsigned i, j, k;
	int unsigned ii, jj, index, jndex;
	double temp;
	
	// Initialize matrices
	for (i = 0; i < SIZE; i++) {
                for (j = 0; j < SIZE; j++) {
                        if (i == j) {
                                A[i][j] = 1.;
                                B[i][j] = 2.;
				C[i][k] = 0.;
                        } else {
                                A[i][j] = 0.;
                                B[i][j] = 2.;
				C[i][j] = 0.;
                        }
                }
        }


	// Matrix mulitply using blocks and multithreading via OpenMP
	// Here, setting 6 threads (can choose)
	#pragma omp parallel for shared(A,B,C) private(i, j, ii, jj, k) schedule(auto) collapse(2) num_threads(6)
        for (i = 0; i < SIZE; i += BLOCK_SIZE) {
                for (j = 0; j < SIZE; j += BLOCK_SIZE) {	
			for (ii = 0; ii < BLOCK_SIZE; ii += 1) {
                                for (jj = 0; jj < BLOCK_SIZE; jj += 1) {
                                        index = ii + i;
                                        jndex = jj + j;
                                        if (index < SIZE && jndex < SIZE) {
						temp = 0.;
						for (k = 0; k < SIZE; k++) {
                                                        temp += A[index][k] * B[k][jndex];
                                                }
						C[index][jndex] = temp;
					}
				}
			}
		}
	}


	// Print product matrix C to check for smaller matrix sizes
/*
	for (i = 0; i < SIZE; i++) {
    		for (j = 0; j < SIZE; j++) {
      			printf ("%.f ", C[i][j]);
    		}
    		printf("\n");
  	}
*/	
	

	// Check number of threads used 
	printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
	printf ( "  Number of threads              = %d\n", omp_get_max_threads ( ) );

        return 0;
}


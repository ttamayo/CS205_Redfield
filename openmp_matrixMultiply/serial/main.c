// Matrix-matrix multiplication
// Serial, naive 3-loop implementation


#include <stdio.h>
#include <stdlib.h>

#define SIZE 8192

int main(void) {

	double A[SIZE][SIZE];
	double B[SIZE][SIZE];
	double C[SIZE][SIZE];

 	int unsigned i, j, k;

	// Initialize matrices
	for (i = 0; i < SIZE; i++) {
		for (j = 0; j < SIZE; j++) {
			if (i == j) {
				A[i][j] = 1.;
        			B[i][j] = 2.;
        			C[i][j] = 0.;
      			} else {
        			A[i][j] = 0.;
        			B[i][j] = 2.;
        			C[i][j] = 0.;
      			}
    		}
  	}


	// Multiply: C = A * B 
  	for (i = 0; i < SIZE; i++) {
    		for (j = 0; j < SIZE; j++) {
      			for (k = 0; k < SIZE; k++) {
        			C[i][j] += A[i][k] * B[k][j];
      			}
    		}
  	}


	// Print C
/*
  for (i = 0; i < SIZE; i++) {
    for (j = 0; j < SIZE; j++) {
      printf ("%.f ", C[i][j]);
    }
    printf("\n");
  }  
*/

  return 0;
}




// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 4
#define dt 1.0
#define number_of_steps 10

/**********************************************************************/


int main(void) {
	printf("# starting\n");
	int unsigned i, j, k;
	double tic, toc;
	int SIZE;
	SIZE = NSITES + 2;
	printf("# SIZE %d\n", SIZE);



	double *rho[SIZE];
	for (i = 0; i < SIZE; i++) {
		rho[i] = (double *) malloc(SIZE * sizeof(double));
	}


	#pragma acc kernels
	for (i = 0; i < SIZE; i++) {
		for (j = 0; j < SIZE; j++) {
			rho[i][j] = 1.;
		}
	}


	double *test;
	test = (double *) malloc(sizeof(double) * SIZE*SIZE);

	#pragma acc kernels 
	#pragma acc loop independent
	for (i = 0; i < SIZE; i++) {
		#pragma acc loop independent
		for (j = 0; j < SIZE; j++) {
			test[j + i * SIZE] = 2.;
		}
	}



/* --------------------------------------------------------- */

	for (i = 0; i < SIZE; i++) {
		for (j = 0; j < SIZE; j++) {
			printf("%.5f  ", rho[i][j]);
		}
		printf("\n");
	}

	printf("\n\n");


        for (i = 0; i < SIZE; i++) {
                for (j = 0; j < SIZE; j++) {
                        printf("%.5f  ", test[i + j * SIZE]);
                }
                printf("\n");
        }

	return 0;
}

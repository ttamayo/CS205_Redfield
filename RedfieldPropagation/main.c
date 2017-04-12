

// external modules
#include <stdio.h>
#include <stdlib.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define SIZE 4

/**********************************************************************/


int main(void) {
	int unsigned i, j, k;

	double *A;
	A = (double *) malloc(sizeof(double) * SIZE * SIZE);

	gen_test_hamiltonian(A);
	print_matrix_real(A, SIZE);
	printf("--\n");

	double *D;
	D = (double *) malloc(sizeof(double) * SIZE);

	printf("set up matrices\n");

	// make A the eigenvector matrix and D the eigenvalue array
	diagonalize(A, D, SIZE);

	for (i = 0; i < SIZE; i++) {
		printf("%.1f ", D[i]);
	}

	printf("diagonalized matrix\n");

	double *gammas;
	gammas = (double *) malloc(sizeof(double) * SIZE * SIZE * SIZE);

	double *params;
	params = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	params[0] = 50;
	params[1] = 35;
	params[2] = 0;
	params[3] = 50;
	params[4] = 35;
	params[5] = 0;
	double *links_to_loss;
	links_to_loss = (double *) malloc(sizeof(double) * SIZE);
	for (i = 0; i < SIZE; i++) 
		links_to_loss[1] = 0.0005;
	double *links_to_target;
	links_to_target = (double *) malloc(sizeof(double) * SIZE);
	links_to_target[1] = 0.005;

	printf("set up rate calculation\n");

	get_rates(gammas, params, D, 2, 4);


	double *V;
	V = (double *) malloc(sizeof(double) * SIZE * SIZE);
	printf("getting V\n");
	get_V(V, A, 1, 1, SIZE);
	printf("printing V\n");
	print_matrix_real(V, SIZE);


	printf("--\n");
	for (i = 0; i < SIZE; i++) {
		for (j = 0; j < SIZE; j++) {
			for (k = 0; k < SIZE; k++) {
				printf("%.5f ", gammas[i + j * SIZE + k * SIZE * SIZE]);
			}
			printf("\n");
		}
		printf("--\n");
	}

	printf("calculated rates\n");

	return 0;

	// FIXME
	// not sure why this should be after the return command
	// putting this line before the return raises an error 
//	free((void*) A);

}
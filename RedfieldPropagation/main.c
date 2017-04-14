

// external modules
#include <stdio.h>
#include <stdlib.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define SIZE 4
#define dt 0.01

/**********************************************************************/


int main(void) {
	int unsigned i, j, k;

	double *A, *D, *V, *gammas, *params, *links_to_loss, *links_to_target;
	A               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	D               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	V               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gammas          = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE);
	params          = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);


	gen_test_hamiltonian(A);
	gen_test_spec_densities(params);
	gen_test_links(links_to_loss, links_to_target, SIZE);

	double *rho_real, *rho_imag, *helper_matrix;
	rho_real      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	rho_imag      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper_matrix = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gen_zero_matrix_complex(rho_real, rho_imag, SIZE);
	rho_real[1 + 1 * SIZE] = 1.;


	// diagonalize hamiltonian
	diagonalize(A, D, SIZE);
	// get rates in site basis
	get_rates(gammas, params, D, 2, 4);
	// rotate rho into exciton basis
	rotate(rho_real, A, SIZE);
	rotate(rho_imag, A, SIZE);

// this is the first propagation step - manually at this point

	int unsigned step, number_of_steps;
	number_of_steps = 50000;

	for (step = 0; step < number_of_steps; step++){

		// get the commutator
		double *comm_real, *comm_imag;
		comm_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
		comm_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);

//		print_matrix_real(rho_real, SIZE);
//		print_matrix_real(rho_imag, SIZE);

		hamiltonian_commutator(rho_real, rho_imag, D, comm_real, comm_imag, SIZE);

		double *lindblad_real, *lindblad_imag;
		lindblad_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
		lindblad_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
		gen_zero_matrix_complex(lindblad_real, lindblad_imag, SIZE);
		lindblad_operator(rho_real, rho_imag, gammas, A, lindblad_real, lindblad_imag, SIZE);

		matrix_mul_scalar(comm_real, dt, SIZE);
		matrix_mul_scalar(comm_imag, dt, SIZE);
		matrix_mul_scalar(lindblad_real, dt, SIZE);
		matrix_mul_scalar(lindblad_imag, dt, SIZE);

//		printf("--- first part ---\n");
//		print_matrix_real(comm_real, SIZE);
//		print_matrix_real(comm_imag, SIZE);
//		printf("------------------\n");


//		printf("--- second part ---\n");
//		print_matrix_real(lindblad_real, SIZE);
//		print_matrix_real(lindblad_imag, SIZE);
//		printf("------------------\n");

		matrix_add_real(rho_real, comm_real, rho_real, SIZE);
		matrix_add_real(rho_real, lindblad_real, rho_real, SIZE);
		matrix_add_real(rho_imag, comm_imag, rho_imag, SIZE);
		matrix_add_real(rho_imag, lindblad_imag, rho_imag, SIZE);

		transpose(A, SIZE);
		rotate(rho_real, A, SIZE);
		rotate(rho_imag, A, SIZE);
		transpose(A, SIZE);

//		printf("new rho:\n");
//		print_matrix_real(rho_real, SIZE);
//		printf("... and ...\n");
//		print_matrix_real(rho_imag, SIZE);
		printf("%d %.10f %.10f\n", step, rho_real[1 + 1 * SIZE], rho_real[2 + 2 * SIZE]);
		rotate(rho_real, A, SIZE);
		rotate(rho_imag, A, SIZE);
	
	}

//	print_matrix_real(comm_real, SIZE);
//	printf("===\n");
//	print_matrix_real(comm_imag, SIZE);
//	printf("===\n");


//	get_V(V, A, 1, 1, SIZE);
//	printf("printing V\n");
//	print_matrix_real(V, SIZE);


//	printf("--\n");
//	for (i = 0; i < SIZE; i++) {
//		for (j = 0; j < SIZE; j++) {
//			for (k = 0; k < SIZE; k++) {
//				printf("%.5f ", gammas[i + j * SIZE + k * SIZE * SIZE]);
//			}
//			printf("\n");
//		}
//		printf("--\n");
//	}
//
//	printf("calculated rates\n");

	return 0;

	// FIXME
	// not sure why this should be after the return command
	// putting this line before the return raises an error 
//	free((void*) A);

}
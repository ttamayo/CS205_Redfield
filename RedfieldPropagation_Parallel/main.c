

// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 8
#define dt 0.025

/**********************************************************************/


int main(void) {
	printf("starting\n");
	int unsigned i, j, k;
	double tic, toc;
	int SIZE;
	SIZE = NSITES + 2;

//	exit(1);

	double *A, *D, *V, *gammas, *params, *links_to_loss, *links_to_target;
	A               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	D               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	V               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gammas          = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE);
	params          = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);

//	exit(1);

	gen_random_hamiltonian_real(A, SIZE);
	gen_test_spec_densities(params, NSITES);
	gen_test_links(links_to_loss, links_to_target, SIZE);

	double *rho_real, *rho_imag, *helper_matrix;
	rho_real      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	rho_imag      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper_matrix = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gen_zero_matrix_complex(rho_real, rho_imag, SIZE);
	rho_real[1 + 1 * SIZE] = 1.;


//	print_matrix_real(rho_real, SIZE);

	// diagonalize hamiltonian
	diagonalize(A, D, SIZE);
	// get rates in site basis
	get_rates(gammas, params, D, NSITES, SIZE);
	// rotate rho into exciton basis
	rotate(rho_real, A, SIZE);
	rotate(rho_imag, A, SIZE);

// this is the first propagation step - manually at this point

	int unsigned step, number_of_steps;
	number_of_steps = 1000;

	tic = clock();

//	#pragma acc kernels loop gang, vector
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
		lindblad_operator(rho_real, rho_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, SIZE);


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

//		printf("%d ", step);
//		for (i = 0; i < SIZE; i++) {
//			printf("%.10f ", rho_real[i + i * SIZE]);
//		}
//		printf("\n");

		rotate(rho_real, A, SIZE);
		rotate(rho_imag, A, SIZE);
	
		free((void*) comm_real);
		free((void*) comm_imag);
		free((void*) lindblad_real);
		free((void*) lindblad_imag);
//
	}
//
	toc = clock();

    // Analyze time elapsed
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n# RESULTS:\n");
    printf("# --------\n");
    printf("# Time Elapsed: %f seconds\n\n", time_spent);

    return 0;

	// FIXME
	// not sure why this should be after the return command
	// putting this line before the return raises an error 
//	free((void*) A);

}

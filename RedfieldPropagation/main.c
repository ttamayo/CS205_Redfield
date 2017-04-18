

// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 6
#define dt 1.0

/**********************************************************************/


int main(void) {
	int unsigned i, j, k;
	double tic, toc;
	int SIZE;
	SIZE = NSITES + 2;

	double *A, *D, *V, *gammas, *params, *links_to_loss, *links_to_target;
	A               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	D               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	V               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gammas          = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE);
	params          = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);



	gen_random_hamiltonian_real(A, SIZE);
	gen_test_spec_densities(params, NSITES);
	gen_test_links(links_to_loss, links_to_target, SIZE);

	double *rho_real, *rho_imag, *helper_matrix;
	rho_real      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	rho_imag      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper_matrix = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gen_zero_matrix_complex(rho_real, rho_imag, SIZE);
	rho_real[1 + 1 * SIZE] = 1.;

//	print_matrix_real(A, SIZE);
//	exit(1);

//	print_matrix_real(rho_real, SIZE);

	// diagonalize hamiltonian
	diagonalize(A, D, SIZE);
	// get rates in site basis
	get_rates(gammas, params, D, NSITES, SIZE);
	// rotate rho into exciton basis
	rotate(rho_real, A, SIZE);
	rotate(rho_imag, A, SIZE);


	double *comm_real, *comm_imag;
	comm_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	comm_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	double *lindblad_real, *lindblad_imag;
	lindblad_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	lindblad_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);

	double *k1_real, *k2_real, *k3_real, *k4_real, *h1_real;
	double *k1_imag, *k2_imag, *k3_imag, *k4_imag, *h1_imag;
	k1_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k1_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k2_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k2_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k3_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k3_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k4_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	k4_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	h1_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	h1_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);


	int unsigned step, number_of_steps;
	number_of_steps = 1000;

	tic = clock();

	for (step = 0; step < number_of_steps; step++){

		// implementing a 4th order runge kutta scheme here
		// get k1
		get_density_update(rho_real, rho_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, SIZE);

//		gen_zero_matrix_complex(lindblad_real, lindblad_imag, SIZE);
		matrix_add_complex(comm_real, comm_imag, lindblad_real, lindblad_imag, k1_real, k1_imag, SIZE);

		// get k2
		matrix_mul_scalar(k1_real, dt/2., SIZE);
		matrix_mul_scalar(k1_imag, dt/2., SIZE);
		matrix_add_complex(rho_real, rho_imag, k1_real, k1_imag, h1_real, h1_imag, SIZE);
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, SIZE);

//		gen_zero_matrix_complex(lindblad_real, lindblad_imag, SIZE);
		matrix_add_complex(comm_real, comm_imag, lindblad_real, lindblad_imag, k2_real, k2_imag, SIZE);

		// get k3
		matrix_mul_scalar(k2_real, dt/2., SIZE);
		matrix_mul_scalar(k2_imag, dt/2., SIZE);
		matrix_add_complex(rho_real, rho_imag, k2_real, k2_imag, h1_real, h1_imag, SIZE);
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, SIZE);

//		gen_zero_matrix_complex(lindblad_real, lindblad_imag, SIZE);
		matrix_add_complex(comm_real, comm_imag, lindblad_real, lindblad_imag, k3_real, k3_imag, SIZE);
		
		// get k4
		matrix_mul_scalar(k3_real, dt, SIZE);
		matrix_mul_scalar(k3_imag, dt, SIZE);
		matrix_add_complex(rho_real, rho_imag, k3_real, k3_imag, h1_real, h1_imag, SIZE);
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, SIZE);

//		gen_zero_matrix_complex(lindblad_real, lindblad_imag, SIZE);
		matrix_add_complex(comm_real, comm_imag, lindblad_real, lindblad_imag, k4_real, k4_imag, SIZE);
		
		// summary:
		// we computed k1 * dt / 2., k2 * dt / 2., k3 * dt, k4

		// now we update the density
		// the update comprises of k1 * dt / 6., k2 * dt / 3., k3 * dt / 3., k4 * dt / 6.
		matrix_mul_scalar(k1_real, 1/3., SIZE);
		matrix_mul_scalar(k1_imag, 1/3., SIZE);
		matrix_mul_scalar(k2_real, 2/3., SIZE);
		matrix_mul_scalar(k2_imag, 2/3., SIZE);
		matrix_mul_scalar(k3_real, 1/3., SIZE);
		matrix_mul_scalar(k3_imag, 1/3., SIZE);
		matrix_mul_scalar(k4_real, dt / 6., SIZE);
		matrix_mul_scalar(k4_imag, dt / 6., SIZE);

		matrix_add_complex(rho_real, rho_imag, k1_real, k1_imag, rho_real, rho_imag, SIZE);
		matrix_add_complex(rho_real, rho_imag, k2_real, k2_imag, rho_real, rho_imag, SIZE);
		matrix_add_complex(rho_real, rho_imag, k3_real, k3_imag, rho_real, rho_imag, SIZE);
		matrix_add_complex(rho_real, rho_imag, k4_real, k4_imag, rho_real, rho_imag, SIZE);

		// done with Runge Kutta step


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
	

	}

	free((void*) comm_real);
	free((void*) comm_imag);
	free((void*) lindblad_real);
	free((void*) lindblad_imag);

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

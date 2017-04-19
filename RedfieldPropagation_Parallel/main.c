

// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 22
#define dt 0.5

/**********************************************************************/


int main(void) {
	printf("# starting\n");
	int unsigned i, j, k;
	double tic, toc;
	int SIZE;
	SIZE = NSITES + 2;

//	exit(1);

	double *A, *D, *gammas, *params, *links_to_loss, *links_to_target;
	A               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	D               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gammas          = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE);
	params          = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);

//	exit(1);

	gen_random_hamiltonian_real(A, SIZE);
	gen_test_spec_densities(params, NSITES);
	gen_test_links(links_to_loss, links_to_target, SIZE);


//	print_matrix_real(A, SIZE);
//	exit(1);

	double *rho_real, *rho_imag, *helper_matrix;
	rho_real      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	rho_imag      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper_matrix = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gen_zero_matrix_complex(rho_real, rho_imag, SIZE);
	rho_real[1 + 1 * SIZE] = 1.;

	double *all_Vs;
	all_Vs = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE*3);


//	print_matrix_real(rho_real, SIZE);

	// diagonalize hamiltonian
	diagonalize(A, D, SIZE);
	// get rates in site basis
	get_V_matrices(all_Vs, A, SIZE);
	
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

//	int unsigned i, j, k;

	int unsigned step, number_of_steps;
	number_of_steps = 5;

	tic = clock();

//	#pragma acc data copyin(gammas[0:SIZE*SIZE*SIZE], D[0:SIZE], A[0:SIZE*SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE]) 
//	#pragma acc data create(k1_real[0:SIZE*SIZE], k1_imag[0:SIZE*SIZE], k2_real[0:SIZE*SIZE], k2_imag[0:SIZE*SIZE], k3_real[0:SIZE*SIZE], k3_imag[0:SIZE*SIZE], k4_real[0:SIZE*SIZE], k4_imag[0:SIZE*SIZE]) 
//	#pragma acc data create(comm_real[0:SIZE*SIZE], comm_imag[0:SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE], h1_real[0:SIZE*SIZE], h1_imag[0:SIZE*SIZE]) 
	#pragma acc data copy(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3])
	for (step = 0; step < number_of_steps; step++){

		// implementing a 4th order runge kutta scheme here

		//=== get k1 ===//	
		get_density_update(rho_real, rho_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);
//		#pragma acc kernels present_or_copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copyin(comm_real[0:SIZE*SIZE], comm_imag[0:SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE]) copy(k1_real[0:SIZE*SIZE], k1_imag[0:SIZE*SIZE], h1_real[0:SIZE*SIZE], h1_imag[0:SIZE*SIZE])
		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				k1_real[i + j * SIZE] = (comm_real[i + j * SIZE] + lindblad_real[i + j * SIZE]) * dt / 2.;
				k1_imag[i + j * SIZE] = (comm_imag[i + j * SIZE] + lindblad_imag[i + j * SIZE]) * dt / 2.;
				h1_real[i + j * SIZE] = rho_real[i + j * SIZE] + k1_real[i + j * SIZE];
				h1_imag[i + j * SIZE] = rho_imag[i + j * SIZE] + k1_imag[i + j * SIZE];
			}
		}

		//=== get k2 ===//
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);
//		#pragma acc kernels present_or_copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copyin(comm_real[0:SIZE*SIZE], comm_imag[0:SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE]) copy(k2_real[0:SIZE*SIZE], k2_imag[0:SIZE*SIZE], h1_real[0:SIZE*SIZE], h1_imag[0:SIZE*SIZE])
		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				k2_real[i + j * SIZE] = (comm_real[i + j * SIZE] + lindblad_real[i + j * SIZE]) * dt / 2.;
				k2_imag[i + j * SIZE] = (comm_imag[i + j * SIZE] + lindblad_imag[i + j * SIZE]) * dt / 2.;
				h1_real[i + j * SIZE] = rho_real[i + j * SIZE] + k2_real[i + j * SIZE];
				h1_imag[i + j * SIZE] = rho_imag[i + j * SIZE] + k2_imag[i + j * SIZE];
			}
		}

		//=== get k3 ===//
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);
//		#pragma acc kernels present_or_copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copyin(comm_real[0:SIZE*SIZE], comm_imag[0:SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE]) copy(k3_real[0:SIZE*SIZE], k3_imag[0:SIZE*SIZE], h1_real[0:SIZE*SIZE], h1_imag[0:SIZE*SIZE])
                #pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				k3_real[i + j * SIZE] = (comm_real[i + j * SIZE] + lindblad_real[i + j * SIZE]) * dt;
				k3_imag[i + j * SIZE] = (comm_imag[i + j * SIZE] + lindblad_imag[i + j * SIZE]) * dt;
				h1_real[i + j * SIZE] = rho_real[i + j * SIZE] + k3_real[i + j * SIZE];
				h1_imag[i + j * SIZE] = rho_imag[i + j * SIZE] + k3_imag[i + j * SIZE];
			}
		}

		//=== get k4 ===//
		get_density_update(h1_real, h1_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);

		//=== combine everyting and update rho ===//
//		#pragma acc kernels copyout(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copyin_or_present(k1_real[0:SIZE*SIZE], k1_imag[0:SIZE*SIZE], k2_real[0:SIZE*SIZE], k2_imag[0:SIZE*SIZE], k3_real[0:SIZE*SIZE], k3_imag[0:SIZE*SIZE])  copyin(comm_real[0:SIZE*SIZE], comm_imag[0:SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				rho_real[i + j * SIZE] += k1_real[i + j * SIZE] / 3. + 2 * k2_real[i + j * SIZE] / 3. + k3_real[i + j * SIZE] / 3. + (comm_real[i + j * SIZE] + lindblad_real[i + j * SIZE]) * dt / 6.;
				rho_imag[i + j * SIZE] += k1_imag[i + j * SIZE] / 3. + 2 * k2_imag[i + j * SIZE] / 3. + k3_imag[i + j * SIZE] / 3. + (comm_imag[i + j * SIZE] + lindblad_imag[i + j * SIZE]) * dt / 6.;
			}
		}


		// done with Runge Kutta step

		//FIXME: In principle, we can keep integrating and rotate back on a different device?

		transpose(A, SIZE);
		rotate(rho_real, A, SIZE);
		rotate(rho_imag, A, SIZE);
		transpose(A, SIZE);

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
}

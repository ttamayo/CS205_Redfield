

// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 48
#define dt 1.0
#define number_of_steps 10

/**********************************************************************/



int main(void) {
	printf("# starting the propagation ...\n");
 
	int unsigned i, j, k, SIZE;
	double tic, toc;

	SIZE = NSITES + 2;
	printf("# ... propagating a %d x %d matrix\n\n", SIZE, SIZE);

	/* allocate matrices */
	// one dimensional
	double *hamiltonian, *links_to_loss, *links_to_target;
	hamiltonian     = (double *) malloc(sizeof(double) * SIZE);
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);
	// two dimensional, squared
	double *eigVects[SIZE], *rho_real[SIZE], *rho_imag[SIZE];
	for (i = 0; i < SIZE; i++){
		eigVects[i] = (double *) malloc(sizeof(double) * SIZE);
		rho_real[i] = (double *) malloc(sizeof(double) * SIZE);
		rho_imag[i] = (double *) malloc(sizeof(double) * SIZE);
	}
	// two dimensional, rectangular
	double *params[SIZE - 2];
	for (i = 0; i < SIZE - 2; i++) {
		params[i] = (double *) malloc(sizeof(double) * 3);
	}
	// three dimensional
	double **gammas[SIZE];
	for (i = 0; i < SIZE; i++) {
		gammas[i] = (double **) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			gammas[i][j] = (double *) malloc(sizeof(double) * SIZE);
		}
	}	
	// four dimensional, rectangular
	double ***all_Vs[SIZE];
	for (i = 0; i < SIZE; i++) {
		all_Vs[i] = (double ***) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			all_Vs[i][j] = (double **) malloc(sizeof(double) * SIZE);
			for (k = 0; k < SIZE; k++) {
				all_Vs[i][j][k] = (double *) malloc(sizeof(double) * 3);
			}
		}
	}

	/* initialize matrices */
	gen_random_hamiltonian_real(eigVects, SIZE);
	gen_test_spec_densities(params, NSITES);
	gen_test_links(links_to_loss, links_to_target, SIZE);
	gen_zero_matrix(rho_real, SIZE);
	gen_zero_matrix(rho_imag, SIZE);
	rho_real[1][1] = 1.;

	/* preprocess matrices */
	diagonalize(eigVects, hamiltonian, SIZE);
	get_rates(gammas, params, hamiltonian, SIZE);		// produces a bunch of nans
	get_V_matrices(all_Vs, eigVects, SIZE);

	rotate(rho_real, eigVects, SIZE);
	rotate(rho_imag, eigVects, SIZE);

	/* allocate more matrices */
	// two dimensional, squared
	double *comm_real[SIZE], *comm_imag[SIZE], *lindblad_real[SIZE], *lindblad_imag[SIZE];
	for (i = 0; i < SIZE; i++){
		comm_real[i] = (double *) malloc(sizeof(double) * SIZE);
		comm_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		lindblad_real[i] = (double *) malloc(sizeof(double) * SIZE);
		lindblad_imag[i] = (double *) malloc(sizeof(double) * SIZE);
	}
	double *k1_real[SIZE], *k2_real[SIZE], *k3_real[SIZE], *k4_real[SIZE], *h1_real[SIZE];
	double *k1_imag[SIZE], *k2_imag[SIZE], *k3_imag[SIZE], *k4_imag[SIZE], *h1_imag[SIZE];
	for (i = 0; i < SIZE; i++) {
		k1_real[i] = (double *) malloc(sizeof(double) * SIZE);
		k2_real[i] = (double *) malloc(sizeof(double) * SIZE);
		k3_real[i] = (double *) malloc(sizeof(double) * SIZE);
		k4_real[i] = (double *) malloc(sizeof(double) * SIZE);
		h1_real[i] = (double *) malloc(sizeof(double) * SIZE);
		k1_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		k2_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		k3_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		k4_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		h1_imag[i] = (double *) malloc(sizeof(double) * SIZE);
	} 


	// start propagation //

	int unsigned step;
	tic = clock();

//	#pragma acc data copyin(hamiltonian[0:SIZE])   copy(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])  //  create(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], k4_real[0:SIZE][0:SIZE], k4_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
	#pragma acc data copyin(hamiltonian[0:SIZE], links_to_target[0:SIZE], links_to_loss[0:SIZE], gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][3])
	for (step = 0; step < number_of_steps; step++) {

		// propagate in a 4th order runge kutta scheme

		//=== get k1 ===//
		
		get_density_update(rho_real, rho_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);

		#pragma acc kernels
		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				k1_real[i][j] = (comm_real[i][j] + lindblad_real[i][j]) * dt / 2.;
				k1_imag[i][j] = (comm_imag[i][j] + lindblad_imag[i][j]) * dt / 2.;
				h1_real[i][j] = rho_real[i][j] + k1_real[i][j];
				h1_imag[i][j] = rho_imag[i][j] + k1_imag[i][j];

			}
		}

		//=== get k2 ===//

//		printf("getting k2\n");
		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);

		#pragma acc kernels
		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				k2_real[i][j] = (comm_real[i][j] + lindblad_real[i][j]) * dt / 2.;
				k2_imag[i][j] = (comm_imag[i][j] + lindblad_imag[i][j]) * dt / 2.;
				h1_real[i][j] = rho_real[i][j] + k2_real[i][j];
				h1_imag[i][j] = rho_imag[i][j] + k2_imag[i][j];

			}
		}


		//=== get k3 ===//

//		printf("getting k3\n");
		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);

		#pragma acc kernels
		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				k3_real[i][j] = (comm_real[i][j] + lindblad_real[i][j]) * dt ;
				k3_imag[i][j] = (comm_imag[i][j] + lindblad_imag[i][j]) * dt ;
				h1_real[i][j] = rho_real[i][j] + k3_real[i][j];
				h1_imag[i][j] = rho_imag[i][j] + k3_imag[i][j];

			}
		}	

		//=== get k4 ===//

//		printf("getting k4\n");
		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);

		#pragma acc kernels
		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				rho_real[i][j] += k1_real[i][j] / 3. + 2 * k2_real[i][j] / 3. + k3_real[i][j] / 3. + (comm_real[i][j] + lindblad_real[i][j]) * dt / 6.;
				rho_imag[i][j] += k1_imag[i][j] / 3. + 2 * k2_imag[i][j] / 3. + k3_imag[i][j] / 3. + (comm_imag[i][j] + lindblad_imag[i][j]) * dt / 6.;
			}
		}

		transpose(eigVects, SIZE);
		rotate(rho_real, eigVects, SIZE);
		rotate(rho_imag, eigVects, SIZE);
		transpose(eigVects, SIZE);


		printf("%d ", step);
        for (i = 0; i < SIZE; i++) {
            printf("%.10f ", rho_real[i][i]);
        }
        printf("\n");

        rotate(rho_real, eigVects, SIZE);
        rotate(rho_imag, eigVects, SIZE);


	}

// 		free((void*) comm_real);
//        free((void*) comm_imag);
//        free((void*) lindblad_real);
//        free((void*) lindblad_imag);

        toc = clock();

    // Analyze time elapsed
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\n# RESULTS:\n");
    printf("# --------\n");
    printf("# Time Elapsed: %f seconds\n\n", time_spent);



	return 0;
}

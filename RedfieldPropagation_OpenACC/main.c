

// external modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define BLOCK_SIZE 16
#define NSITES 62
#define dt 1.0
#define number_of_steps 10

/**********************************************************************/



int main(void) {
	printf("# starting the propagation ...\n");
 
	int unsigned i, j, k, l, SIZE;
	double tic, toc, total;

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
	// three dimensional, squared
	double **reduction_intermediates[SIZE];
	for (i = 0; i < SIZE; i++) {
		reduction_intermediates[i] = (double **) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			reduction_intermediates[i][j] = (double *) malloc(sizeof(double) * SIZE);
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

	// even more helpers
//	double *V[SIZE];
//	for (i = 0; i < SIZE; i++) {
//		V[i] = (double *) malloc(sizeof(double) * SIZE);
//	}
//	double *first_real[SIZE], *second_real[SIZE], *helper[SIZE];
//	for (i = 0; i < SIZE; i++) {
//		first_real[i] = (double *) malloc(sizeof(double) * SIZE);
//		second_real[i] = (double *) malloc(sizeof(double) * SIZE);
//		helper[i] = (double *) malloc(sizeof(double) * SIZE);
//	}

	double **V[SIZE], **first[SIZE], **second[SIZE], **helper[SIZE];
	for (i = 0; i < SIZE; i++) {
		V[i]      = (double **) malloc(sizeof(double) * SIZE);
		first[i]  = (double **) malloc(sizeof(double) * SIZE);
		second[i] = (double **) malloc(sizeof(double) * SIZE);
		helper[i] = (double **) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			V[i][j]      = (double *) malloc(sizeof(double) * SIZE);
			first[i][j]  = (double *) malloc(sizeof(double) * SIZE);
			second[i][j] = (double *) malloc(sizeof(double) * SIZE);
			helper[i][j] = (double *) malloc(sizeof(double) * SIZE);
		}
	}


	// start propagation //
	int unsigned ii, jj;
	int unsigned index, jndex, kndex;
	int unsigned step;
	tic = clock();

	printf("started propagation\n");

//	#pragma acc data copyin(hamiltonian[0:SIZE], eigVects[0:SIZE][0:SIZE], links_to_target[0:SIZE], links_to_loss[0:SIZE], gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3])   copyin(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE])   copyin(h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE], k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE])    copyin(V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE])     copy(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])    create(reduction_intermediates[0:SIZE][0:SIZE][0:SIZE][0:SIZE][0:SIZE])
//	#pragma acc data copyin(hamiltonian[0:SIZE])

	#pragma acc data copyin(hamiltonian[0:SIZE], gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], links_to_target[0:SIZE], links_to_loss[0:SIZE])  
	#pragma acc data copyin(V[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc data copyin(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
	#pragma acc data copyin(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE])
	#pragma acc data create(reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc data copyin(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])
	for (step = 0; step < number_of_steps; step++) {

		// propagate in a 4th order runge kutta scheme

		//=== get k1 ===//
	
//		#pragma acc kernels	
//		{
		get_density_update(rho_real, rho_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);

		#pragma acc kernels     present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
		#pragma acc loop independent collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
		                       	         	k1_real[index][jndex] = (comm_real[index][jndex] + lindblad_real[index][jndex]) * dt / 2.;
       		                       			k1_imag[index][jndex] = (comm_imag[index][jndex] + lindblad_imag[index][jndex]) * dt / 2.;
                		                	h1_real[index][jndex] = rho_real[index][jndex] + k1_real[index][jndex];
                        		        	h1_imag[index][jndex] = rho_imag[index][jndex] + k1_imag[index][jndex];
						}
					}					
				}
			}		

		//=== get k2 ===//
		
//		printf("getting k2\n");
//		#pragma acc update host(h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])

		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);
                
		#pragma acc kernels present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
		#pragma acc loop independent collapse(2)
                for (ii = 0; ii < SIZE; ii += BLOCK_SIZE)
                        for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
                                #pragma acc loop independent collapse(2)
                                for (i = 0; i < BLOCK_SIZE; i++) {
                                        for (j = 0; j < BLOCK_SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        k2_real[index][jndex] = (comm_real[index][jndex] + lindblad_real[index][jndex]) * dt / 2.;
                                                        k2_imag[index][jndex] = (comm_imag[index][jndex] + lindblad_imag[index][jndex]) * dt / 2.;
                                                        h1_real[index][jndex] = rho_real[index][jndex] + k2_real[index][jndex];
                                                        h1_imag[index][jndex] = rho_imag[index][jndex] + k2_imag[index][jndex];
                                                }
                                        } 
                                }
                        }

		//=== get k3 ===//

		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);

		#pragma acc kernels present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
                #pragma acc loop independent collapse(2)
                for (ii = 0; ii < SIZE; ii += BLOCK_SIZE)
                        for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
                                #pragma acc loop independent collapse(2)
                                for (i = 0; i < BLOCK_SIZE; i++) {
                                        for (j = 0; j < BLOCK_SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        k3_real[index][jndex] = (comm_real[index][jndex] + lindblad_real[index][jndex]) * dt;
                                                        k3_imag[index][jndex] = (comm_imag[index][jndex] + lindblad_imag[index][jndex]) * dt;
                                                        h1_real[index][jndex] = rho_real[index][jndex] + k3_real[index][jndex];
                                                        h1_imag[index][jndex] = rho_imag[index][jndex] + k3_imag[index][jndex];
                                                }
                                        }
                                }
                        }

		//=== get k4 ===//

		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);

                #pragma acc kernels present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE], k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE])
                #pragma acc loop independent collapse(2)
                for (ii = 0; ii < SIZE; ii += BLOCK_SIZE)
                        for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
                                #pragma acc loop independent collapse(2)
                                for (i = 0; i < BLOCK_SIZE; i++) {
                                        for (j = 0; j < BLOCK_SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
							rho_real[index][jndex] += k1_real[index][jndex] / 3. + 2 * k2_real[index][jndex] / 3. + k3_real[index][jndex] / 3. + (comm_real[index][jndex] + lindblad_real[index][jndex]) * dt / 6.;
							rho_imag[index][jndex] += k1_imag[index][jndex] / 3. + 2 * k2_imag[index][jndex] / 3. + k3_imag[index][jndex] / 3. + (comm_imag[index][jndex] + lindblad_imag[index][jndex]) * dt / 6.;
						}
					}
				}
			}
//		} // end kernels

		#pragma acc update host(rho_real[0:SIZE][0:SIZE])


		transpose(eigVects, SIZE);
		rotate(rho_real, eigVects, SIZE);
		rotate(rho_imag, eigVects, SIZE);
		transpose(eigVects, SIZE);


	total = 0;
	printf("%d ", step);
        for (i = 0; i < SIZE; i++) {
            printf("%.10f ", rho_real[i][i]);
		total += rho_real[i][i];
        }
	printf("%.10f ", total);
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
// acc data present(h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE], k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], rho_real[0:SIZE], rho_imag[0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE])



	return 0;
}

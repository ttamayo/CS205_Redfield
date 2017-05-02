/*
 * Main program for propagating an excitonic system based on 
 * the secular Redfield approximation
 * 
 * Exciton Hamiltonians of size NSITES x NSITES are randomly drawn
 * 
 * Trap and sink states are explicitly modeled as 
 * two additional states in the exciton Hamiltonian
 *
 * The system is propagated in a 4th order Runge Kutta scheme
 *
 * Exciton populations are recorded at every iteration step
 *
 * author: Florian Hase
 *
 */

// loading all modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "headers.h"

//********************************************************************//
// defining global variables 

#define NSITES 62			 // number of excitonic sites to be modeled
#define dt 1.0				 // integration time step in fs
#define number_of_steps 10   // number of integration time steps
#define BLOCK_SIZE 16        // size for blocked matrix operations

//********************************************************************//
// main routine

int main(void) {
	printf("# starting the propagation ...\n");

	// we start with initializing a number of helper integers and floats
	int unsigned i, j, k, l, SIZE;
	double tic, toc, total;

	// compute the number of states in the Hamiltonian
	// ... remember: sinks and traps are modeled explicitly -> 2 more states
	SIZE = NSITES + 2;
	printf("# ... propagating a %d x %d hamiltonian\n\n", SIZE, SIZE);

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... hamiltonian, eigenvectors, eigenvalues
	// ... rates to sinks and target
	// ... real and imaginary part of the density matrix

	double *hamiltonian, *links_to_loss, *links_to_target;			// one dimensional arrays
	hamiltonian     = (double *) malloc(sizeof(double) * SIZE);
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);

	double *eigVects[SIZE], *rho_real[SIZE], *rho_imag[SIZE];		// two dimensional arrays
	for (i = 0; i < SIZE; i++){
		eigVects[i] = (double *) malloc(sizeof(double) * SIZE);
		rho_real[i] = (double *) malloc(sizeof(double) * SIZE);
		rho_imag[i] = (double *) malloc(sizeof(double) * SIZE);
	}

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... spectral density parameters (describe coupling to environment)
	// ... transition rates between excitonic states

	double *params[SIZE - 2];										// two dimensional array
	for (i = 0; i < SIZE - 2; i++) {
		params[i] = (double *) malloc(sizeof(double) * 3);
	}

	double **gammas[SIZE];											// three dimensional array
	for (i = 0; i < SIZE; i++) {
		gammas[i] = (double **) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			gammas[i][j] = (double *) malloc(sizeof(double) * SIZE);
		}
	}	

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... storing transition matrices V
	// ... storing intermediate results

	double ***all_Vs[SIZE];											// four dimensional array
	for (i = 0; i < SIZE; i++) {
		all_Vs[i] = (double ***) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			all_Vs[i][j] = (double **) malloc(sizeof(double) * SIZE);
			for (k = 0; k < SIZE; k++) {
				all_Vs[i][j][k] = (double *) malloc(sizeof(double) * 3);
			}
		}
	}

	double **reduction_intermediates[SIZE];							// three dimensional array
	for (i = 0; i < SIZE; i++) {
		reduction_intermediates[i] = (double **) malloc(sizeof(double) * SIZE);
		for (j = 0; j < SIZE; j++) {
			reduction_intermediates[i][j] = (double *) malloc(sizeof(double) * SIZE);
		}
	}

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... more intermediate results 
	// ... for computing the density matrix update
	// ... the density matrix update
	// ... the Runge Kutta steps

	// two dimensional arrays
	double *comm_real[SIZE], *comm_imag[SIZE], *lindblad_real[SIZE], *lindblad_imag[SIZE];
	for (i = 0; i < SIZE; i++){
		comm_real[i] = (double *) malloc(sizeof(double) * SIZE);
		comm_imag[i] = (double *) malloc(sizeof(double) * SIZE);
		lindblad_real[i] = (double *) malloc(sizeof(double) * SIZE);
		lindblad_imag[i] = (double *) malloc(sizeof(double) * SIZE);
	}
	// two dimensional arrays
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

	// three dimensional arrays
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

	//------------------------------------------------------------------//

	// now we can set up our system

	// first, generate the hamiltonian
	gen_random_hamiltonian_real(eigVects, SIZE);
	// generate spectral density parameters
	gen_test_spec_densities(params, NSITES);
	// generate rates to sinks and traps
	gen_test_links(links_to_loss, links_to_target, SIZE);
	// generate the initial density matrix - we start the propagation at site 1
	gen_zero_matrix(rho_real, SIZE);
	gen_zero_matrix(rho_imag, SIZE);
	rho_real[1][1] = 1.;

	// diagonalize the hamiltonian
	diagonalize(eigVects, hamiltonian, SIZE);
	// get transition rates between exciton states
	get_rates(gammas, params, hamiltonian, SIZE);
	// get transition matrices 
	get_V_matrices(all_Vs, eigVects, SIZE);

	// rotate the density matrix into the exciton eigenbasis
	rotate(rho_real, eigVects, SIZE);
	rotate(rho_imag, eigVects, SIZE);


	// start propagation //
	int unsigned ii, jj;
	int unsigned index, jndex, kndex;
	int unsigned step;
	tic = clock();

	printf("started propagation\n");

	// copy all data to the GPU
	// ... note: once all the data is on the GPU, only the density matrix needs to be communicated between GPU and CPU
	#pragma acc data copyin(hamiltonian[0:SIZE], gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], links_to_target[0:SIZE], links_to_loss[0:SIZE])  
	#pragma acc data copyin(V[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc data copyin(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
	#pragma acc data copyin(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE])
	#pragma acc data create(reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc data copyin(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])
	// propagate for the number of specified integration steps
	for (step = 0; step < number_of_steps; step++) {
		// a 4th order Runge Kutta scheme is used for propagation
		// as a compromise of computational demand and numerical accuracy

		//=== getting Runge Kutta k1 ===//
		// ... compute density matrix update
		get_density_update(rho_real, rho_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);
		// ... compute k1 in blocked matrix operations
		#pragma acc kernels present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
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

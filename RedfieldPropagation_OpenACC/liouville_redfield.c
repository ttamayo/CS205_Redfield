/* 
 * Calculates the contributions of density matrix updates
 * in the secular Redfield approximation
 * 
 * Computes the commutator -i [H, rho] / hbar
 *
 * and the Lindblad operator 
 * sum_{m,omega} V_m(omega) rho V^dagg_m(omega) - {V^dagg_m(omega) V_m(omega), rho} / 2
 * 
 * in two separate routines and combines the results 
 *
 * computing the Lindblad operator is the most time consuming part of the Redfield propagation
 * 
 * most matrix operations were simplified using properties of the involved matrices
 * ... hamiltonian is real and symmetric, eigenvectors are real and orthogonal, etc.
 *
 * author: Florian Hase
 *
 */


//********************************************************************//
// loading all modules
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "headers.h"


//********************************************************************//
// defining global variables

#define BLOCK_SIZE 8

//********************************************************************//
// defining physical units
// ... energies     ->  cm^-1
// ... time         ->  fs
// ... temperature  ->  K

// unit conversions
#define cm1_to_fs1  1. / 33356.40952
#define fs1_to_cm1  33356.40952
#define eV_to_cm1   8065.54429
#define s_to_fs     1.e15

// physical constants
#define HBAR        6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV    1. / (6.582119514e-16*eV_to_cm1*s_to_fs)
#define KB          8.6173303e-5*eV_to_cm1
#define PI          3.1415926532897932384626433832
#define T           300.0 // FIXME: temperature should be a parameter
#define BETA        1. / (8.6173303e-5*eV_to_cm1 * 300.0)


//********************************************************************//
// computes the commutator between the density matrix and the hamiltonian
// hamiltonian is in diagonal form, which allows to reduce the matrix matrix multiplications
// ... to matrix vector multiplications instead
// parallelized on the gang, worker and vector level
// this step is not performance critical
//   rho_real    ... real part of the density matrix
//   rho_imag    ... imaginary part of the density matrix
//   hamiltonian ... diagonalized hamiltonan
//   comm_real   ... real part of the result
//   comm_imag   ... imaginary part of the result
//   N           ... size of the involved matrices / vectors

void hamiltonian_commutator(double **rho_real, double **rho_imag, double *hamiltonian, double **comm_real, double **comm_imag, int N) {
	int unsigned i, j, ii, jj, index, jndex;
	// blocked matrix vector operations
	// all variables are present on the GPU
	#pragma acc kernels present(hamiltonian[0:N])    present(comm_real[0:N][0:N], comm_imag[0:N][0:N])     present(rho_real[0:N][0:N], rho_imag[0:N][0:N])
	#pragma acc loop independent collapse(2)
    for (ii = 0; ii < N; ii += BLOCK_SIZE)
			for(jj = 0; jj < N; jj += BLOCK_SIZE) {
					#pragma acc loop independent collapse(2)
					for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
									index = ii + i;
									jndex = jj + j;
									if (index < N && jndex < N) {
											comm_real[index][jndex] = (hamiltonian[jndex] - hamiltonian[index]) * rho_imag[index][jndex] * HBAR_INV;
											comm_imag[index][jndex] = (hamiltonian[index] - hamiltonian[jndex]) * rho_real[index][jndex] * HBAR_INV;
									}
							}
					}
			}  
}


//********************************************************************//
// computes the action of the Lindblad operator on the density matrix
// ... in three steps:
// 		 - intra-exciton transitions  (n^6 scaling)
// 		 - loss transitions           (n^4 scaling)
//       - target transitons          (n^4 scaling)
// parallelized on the gang, worker and vector level
// the first step of this function is performance critical
//   rho_real				 ... real part of the density matrix
//   rho_imag				 ... imaginary part of the density matrix
//   gammas					 ... three dimensional tensor containing transition rates
//   eigVects 				 ... eigenvectors of the hamiltonian
//   lindblad_real			 ... real part of the result
//   lindblad_imag			 ... imaginary part of the result
//   links_to_loss			 ... rates to loss state
//   links_to_target		 ... rates to target state
//	 all_Vs					 ... four dimensional tensor containing transition matrices
//   V     					 ... three dimensional tensor for thread block dependent transition matrix storage
//   first, second, helper   ... three dimensional tensors for thread block dependent storage of intermediate results
//   reduction_intermediates ... three dimensional tensor for thread block dependent storage of intermediate results

void lindblad_operator(double **rho_real, double **rho_imag, double ***gammas, double **eigVects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates, int SIZE) {
	// declare helper variables
	int unsigned i, j, k, m, M, N;
	int unsigned ii, jj, kk;
	int unsigned index, jndex, kndex;
	double rate, sum;	
	double iupper, jupper;

	// note: ... for all operations, transpose of matrices are implemented by proper indexing

	// iterate over intra-excitonic transitions (n^6 scaling)
	#pragma acc kernels  num_gangs(128)  num_workers(16)  vector_length(64)  present(gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], V[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE], reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
 	// no parallelization level left for this loop
 	#pragma acc loop
	for (m = 1; m < SIZE - 1; m++) {
		// no parallelization level left for this loop
		#pragma acc loop		
		for (M = 0; M < SIZE; M++) {
			// parallelize over thread blocks
			// ... note: collapsing with outer loops did not improve speed ups for larger matrices
			#pragma acc loop gang independent
			for (N = 0; N < SIZE; N++) {
 
				// extract the transition matrix for this thread block
				// blocked matrix operation parallelized on the worker and vector level
				// using the fast memory of the GPU
				#pragma acc loop  worker  collapse(2) independent
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE)
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)  
						for (i = 0; i < BLOCK_SIZE; i++) {
							for(j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									V[N][index][jndex] = 0;
								}
							}
						}
					}
				V[N][M][N] = all_Vs[M][N][m][1];

			//--> start computing    V rho V^dagg    (first term of Lindblad operator)
				// compute   V rho
				// blocked matrix operation parallelized on the worker and vector level
				// using the fast memory of the GPU
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i; 
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									sum = 0;
									#pragma acc loop
									for (k = 0; k < SIZE; k++) {
										sum += V[N][index][k] * rho_real[k][jndex];
									}
									helper[N][index][jndex] = sum;
								}
							}	
						}
					} 
				
				// now compute   (V rho) V^dagg  with (V rho) stored in helper	
				// blocked matrix operation parallelized on the worker and vector level
				// using the fast memory of the GPU		
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE)
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									sum = 0;
									#pragma acc loop 
									for (k = 0; k < SIZE; k++) {
										sum += helper[N][index][k] * V[N][jndex][k];
									}
									first[N][index][jndex] = sum;
								}
							}
						}
					}

			//--> done with the first term of the Lindblad operator
			//--> on to the next term:   V^dagg V rho / 2.
			//--> ... note: we still have   (V rho)  stored in helper
				// compute V^dagg (V rho) / 2.
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									sum = 0;
									#pragma acc loop seq
									for (k = 0; k < SIZE; k++) {
										sum += V[N][k][index] * helper[N][k][jndex] / 2.;
									}
									second[N][index][jndex] = sum;
								}
							}
						}
					}

			//--> done with the second term of the Lindblad operator
			//--> Let's combine the first two ( V rho V^dagg - V^dagg V rho / 2. ) and store result in first
				// subtract first - second
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector  collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									first[N][index][jndex] = first[N][index][jndex] - second[N][index][jndex];
								}
							}		
						}
					}
	
			//--> now we compute the last term in the Lindblad operator
			//--> which is rho V^dagg V / 2.
				// compute V^dagg V first
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									sum = 0;
									#pragma acc loop 
									for (k = 0; k < SIZE; k++) {
										sum += V[N][k][index] * V[N][k][jndex];
									}
									helper[N][index][jndex] = sum;
								}
							}
						}
					}

				// now compute rho (V^dagg V) / 2. and store result in second
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									sum = 0;
									#pragma acc loop
									for (k = 0; k < SIZE; k++) {			
										sum += rho_real[index][k] * helper[N][k][jndex] / 2.;
									}
									second[N][index][jndex] = sum;
								}
							}
						}
					}

			//--> by now we have also the third term and just need to combine it with the first two
				// calculate first - second
				// store the result in the correct slice of reduction_intermediates
				// for reduction after the operation on multiple thread blocks
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
									reduction_intermediates[index][jndex][N] = gammas[M][N][m] * (first[N][index][jndex] - second[N][index][jndex]);
								}
							}
						}
					}
			
			
			} // end of N loop 			

			// reduce reduction_intermediates and update lindblad_real with intra-exciton transitions
			sum = 0;
			#pragma acc loop independent worker vector collapse(2)			
			for (index = 0; index < SIZE; index++) {
				for (jndex = 0; jndex < SIZE; jndex++) {
					sum = 0;
					#pragma acc loop 
					for (N = 0; N < SIZE; N++) {
						sum += reduction_intermediates[index][jndex][N];
					}
					lindblad_real[index][jndex] += sum;
				}
			}

		} // end of M loop 

	} // end of m loop


	// so far we only computed intra-exciton transitions
	// we still need loss transitions and target transitions
	// these transitions only scale as n^4 and are less critical for performance

	// calculations of these transition follow the same scheme as above
	// except that we don't need to loop over M and N again (which is why we save n^2 here)


	// calculate loss and target transitions in one big kernel
	#pragma acc kernels   num_workers(16) vector_length(64)  present(helper[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], V[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], rho_real[0:SIZE][0:SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], lindblad_real[0:SIZE][0:SIZE])
	{

	//=== getting the losses ===//

	// loop over all states
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		// extract transition matrix
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector, collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							V[m][index][jndex] = all_Vs[index][jndex][m][0];
						}
					}
				}
			}

	//--> start computing    V rho V^dagg    (first term of Lindblad operator)
		// compute   V rho
		// blocked matrix operation parallelized on the worker and vector level
		// using the fast memory of the GPU
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector, collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][index][k] * rho_real[k][jndex];
							}
							helper[m][index][jndex] = sum;
						}
					}
				}
			}

		// now compute   (V rho) V^dagg  with (V rho) stored in helper	
		// blocked matrix operation parallelized on the worker and vector level
		// using the fast memory of the GPU	
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector, collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += helper[m][index][k] * V[m][jndex][k];
							}
							first[m][index][jndex] = sum;
						}
					}
				}
			}


	//--> done with the first term of the Lindblad operator
	//--> on to the next term:   V^dagg V rho / 2.
	//--> ... note: we still have   (V rho)  stored in helper
		// compute V^dagg (V rho) / 2.
		#pragma acc loop independent gang, worker, collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][k][index] * helper[m][k][jndex] / 2.;
							}
							second[m][index][jndex] = sum;
						}
					}
				}
			}


	//--> done with the second term of the Lindblad operator
	//--> Let's combine the first two ( V rho V^dagg - V^dagg V rho / 2. ) and store result in first
		// subtract first - second
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							first[m][index][jndex] = first[m][index][jndex] - second[m][index][jndex];
						}
					}
				}
			}


	//--> now we compute the last term in the Lindblad operator
	//--> which is rho V^dagg V / 2.
		// compute V^dagg V first
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector, collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][k][index] * V[m][k][jndex];
							}
							helper[m][index][jndex] = sum;
						}
					}
				}
			}


		// now compute rho (V^dagg V) / 2. and store result in second
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += rho_real[index][k] * helper[m][k][jndex] / 2.;
							}
							second[m][index][jndex] = sum;
						}
					}
				}
			}

	//--> by now we have also the third term and just need to combine it with the first two
		// calculate first - second
		// multiply by the rate
		// update lindblad_real
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {	
							lindblad_real[index][jndex] += rate * (first[m][index][jndex] - second[m][index][jndex]);
						}
					}
				}
			}

	}  // end of for-loop | not the end of the kernel!


	//=== getting the targets ===//
	// we follow the same protocol again
	// ... note: implementing these steps in separate functions significantly slowed down the code
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];


		// extract transition matrix
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {	
							V[m][index][jndex] = all_Vs[index][jndex][m][2];
						}
					}
				}
			}

	//--> start computing    V rho V^dagg    (first term of Lindblad operator)
		// compute   V rho
		// blocked matrix operation parallelized on the worker and vector level
		// using the fast memory of the GPU
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][index][k] * rho_real[k][jndex];
							}
							helper[m][index][jndex] = sum;
						}
					}
				}
			}

		// now compute   (V rho) V^dagg  with (V rho) stored in helper	
		// blocked matrix operation parallelized on the worker and vector level
		// using the fast memory of the GPU	
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += helper[m][index][k] * V[m][jndex][k];
							}
							first[m][index][jndex] = sum;
						}
					}
				}
			}


	//--> done with the first term of the Lindblad operator
	//--> on to the next term:   V^dagg V rho / 2.
	//--> ... note: we still have   (V rho)  stored in helper
		// compute V^dagg (V rho) / 2.
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][k][index] * helper[m][k][jndex] / 2.;
							}
							second[m][index][jndex] = sum;
						}
					}
				}
			}


	//--> done with the second term of the Lindblad operator
	//--> Let's combine the first two ( V rho V^dagg - V^dagg V rho / 2. ) and store result in first
		// subtract first - second
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {	
							first[m][index][jndex] = first[m][index][jndex] - second[m][index][jndex];
						}
					}
				}
			}



	//--> now we compute the last term in the Lindblad operator
	//--> which is rho V^dagg V / 2.
		// compute V^dagg V first
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += V[m][k][index] * V[m][k][jndex];
							}
							helper[m][index][jndex] = sum;
						}
					}
				}
			}



		// now compute rho (V^dagg V) / 2. and store result in second
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {
							sum = 0;
							#pragma acc loop seq
							for (k = 0; k < SIZE; k++) {
								sum += rho_real[index][k] * helper[m][k][jndex] / 2.;
							}
							second[m][index][jndex] = sum;
						}
					}
				}
			}


	//--> by now we have also the third term and just need to combine it with the first two
		// calculate first - second
		// multiply by the rate
		// update lindblad_real
		#pragma acc loop independent gang, worker collapse(2)
		for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
			for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
				#pragma acc loop independent vector collapse(2)
				for (i = 0; i < BLOCK_SIZE; i++) {
					for (j = 0; j < BLOCK_SIZE; j++) {
						index = ii + i;
						jndex = jj + j;
						if (index < SIZE && jndex < SIZE) {	
							lindblad_real[index][jndex] += rate * (first[m][index][jndex] - second[m][index][jndex]);
						}
					}
				}
			}


	} 
	}

}


//********************************************************************//
// computes the density matrix update in the secular Redfield approximation 
// ... in two steps:
// 		 - Hamiltonian commutator    
//		 - Lindblad operator 
// the second step of this function is performance critical
//   rho_real				 ... real part of the density matrix
//   rho_imag				 ... imaginary part of the density matrix
//   gammas					 ... three dimensional tensor containing transition rates
//   eigVects 				 ... eigenvectors of the hamiltonian
//   comm_real               ... real part of the commutator result
//   comm_imag 				 ... imaginary part of the commutator result
//   lindblad_real			 ... real part of the lindblad result
//   lindblad_imag			 ... imaginary part of the lindblad result
//   links_to_loss			 ... rates to loss state
//   links_to_target		 ... rates to target state
//	 all_Vs					 ... four dimensional tensor containing transition matrices
//   V     					 ... three dimensional tensor for thread block dependent transition matrix storage
//   first, second, helper   ... three dimensional tensors for thread block dependent storage of intermediate results
//   reduction_intermediates ... three dimensional tensor for thread block dependent storage of intermediate results

void get_density_update(double **rho_real, double **rho_imag, double *energies, double **comm_real, double **comm_imag, 
	                    double ***gammas, double **eigvects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates, int N) {

	// get the commutator
	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

	int unsigned i, j;
	// set lindblad result to zero
	#pragma acc kernels present(lindblad_real[0:N][0:N], lindblad_imag[0:N][0:N])
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			lindblad_real[i][j] = 0.;
			lindblad_imag[i][j] = 0.;
		}
	}

	// get lindblad
	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, N);

}

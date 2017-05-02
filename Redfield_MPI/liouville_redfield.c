
#include <stdio.h>
#include <stdlib.h>

#include "headers.h"
//#include "multiple_matrix_operations.h"

/* methods for computing the Liouville operator on the density matrix
 * ... all methods assume operators to be in the exciton basis, unless otherwise specified
 */

#define TILE 16
#define BLOCK_SIZE 8

#define cm1_to_fs1 1. / 33356.40952
#define fs1_to_cm1 33356.40952
#define eV_to_cm1  8065.54429
#define s_to_fs    1.e15

#define HBAR 6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV 1. / (6.582119514e-16*eV_to_cm1*s_to_fs)
#define KB   8.6173303e-5*eV_to_cm1
#define PI   3.1415926532897932384626433832
#define T    300.0 // FIXME: temperature should be a parameter
#define BETA 1. / (8.6173303e-5*eV_to_cm1 * 300.0)

/********************************************************************/

void hamiltonian_commutator(double **rho_real, double **rho_imag, double *hamiltonian, double **comm_real, double **comm_imag, int N) {
	// as dense as possible implementation of the hamiltonian commutator -i * [H, rho] / hbar
	int unsigned i, j, ii, jj, index, jndex;

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



void lindblad_operator(double **rho_real, double **rho_imag, double ***gammas, double **eigVects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates, int SIZE) {
	int unsigned i, j, k, m, M, N;
	int unsigned ii, jj, kk;
	int unsigned index, jndex, kndex;
	double rate, sum;	
	double iupper, jupper;

	#pragma acc kernels  num_gangs(128)  num_workers(16)  vector_length(64)  present_or_copyin(gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], V[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE], reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
  // num_workers(4) vector_length(256)   present(gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE], rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE], reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc loop //independent collapse(2)
	for (m = 1; m < SIZE - 1; m++) {
		#pragma acc loop
		for (M = 0; M < SIZE; M++) {
			#pragma acc loop  gang independent
			for (N = 0; N < SIZE; N++) {
 

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


				//matrix_mul_real(helper, V_dagg, first_real, SIZE);
				// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
			//--------
					
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


				// getting the second part of the equation, V_dagg * V * rho / 2.
			//--------
				//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
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

			//--------

				// subtract first_real - second_real
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
	
			//--------
				//matrix_mul_real(V_dagg, V, helper, SIZE);
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



				// matrix_mul_real(rho, helper, second_real, SIZE) / 2.
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

			//--------
			
				#pragma acc loop independent worker collapse(2)
				for (ii = 0; ii < SIZE; ii += BLOCK_SIZE) 
					for (jj = 0; jj < SIZE; jj += BLOCK_SIZE) {
						#pragma acc loop independent vector collapse(2)
						for (i = 0; i < BLOCK_SIZE; i++) {
							for (j = 0; j < BLOCK_SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {
//									lindblad_real[index][jndex] += gammas[M][N][m] * (first_real[index][jndex] - second_real[index][jndex]);
									reduction_intermediates[index][jndex][N] = gammas[M][N][m] * (first[N][index][jndex] - second[N][index][jndex]);
								}
							}
						}
					}
			
			
			} // end of N loop 			

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


	// big matrix reduction
//	#pragma acc kernels   present(reduction_intermediates[0:SIZE][0:SIZE][0:SIZE][0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE])
//	#pragma acc loop independent worker collapse(2)
//	for (index = 0; index < SIZE; index++) 
//		for (jndex = 0; jndex < SIZE; jndex++) {
//			sum = 0;
//			#pragma acc loop independent vector collapse(2)
//			for (m = 1; m < SIZE - 1; m++) 
//				for (M = 0; M < SIZE; M++)
//					#pragma acc loop seq
//					for (N = 0; N < SIZE; N++) {
//						sum += reduction_intermediates[index][jndex][m][M][N];	
//					}
//			lindblad_real[index][jndex] += sum;
//		}

//	#pragma acc wait(1)


	#pragma acc kernels   num_workers(16)     vector_length(64)        present(helper[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], V[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], rho_real[0:SIZE][0:SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], lindblad_real[0:SIZE][0:SIZE])
	{

	//=== getting the losses ===//
//	#pragma acc loop gang independent 
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		// get V matrix
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


	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
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




	//--------

		// subtract first_real - second_real
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


	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
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



		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
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

	//--------

		// subtract first_real - second_real and multiply by the rate
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
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];


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

		// getting the first part of the equation, V * rho * V_dagg
	//--------
		//matrix_mul_real(V, rho_real, helper, SIZE);
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

		// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
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
	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
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

		//--------

		// subtract first_real - second_real
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



		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
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



		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
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


	//--------

		// subtract first_real - second_real and multiply by the rate
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


	} // end of acc data directive
	}

}

void get_density_update(double **rho_real, double **rho_imag, double *energies, double **comm_real, double **comm_imag, 
	                    double ***gammas, double **eigvects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates, int N) {

//	print_matrix_real(rho_real, N);
	//hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

	int unsigned i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			lindblad_real[i][j] = 0.;
			lindblad_imag[i][j] = 0.;
		}
	}
	return;
	#pragma acc data copyin(lindblad_real[0:N][0:N], lindblad_imag[0:N][0:N])
	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, N);

	#pragma acc update host(lindblad_real[0:N][0:N])
//	print_matrix_real(rho_real, N);
//	exit(1);
//	#pragma acc update device(lindblad_real[0:N][0:N])
//	print_matrix_real(lindblad_real, N);

	return;
}

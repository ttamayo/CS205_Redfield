
#include <stdio.h>
#include <stdlib.h>

#include "headers.h"
//#include "multiple_matrix_operations.h"

/* methods for computing the Liouville operator on the density matrix
 * ... all methods assume operators to be in the exciton basis, unless otherwise specified
 */

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

void hamiltonian_commutator(double * rho_real, double * rho_imag, double * hamiltonian, double * comm_real, double * comm_imag, int N) {
	// as dense as possible implementation of the hamiltonian commutator -i * [H, rho] / hbar
	int unsigned i, j;

	#pragma acc loop independent
	for (i = 0; i < N; i++) {
		#pragma acc loop independent private(comm_real[0:N*N]) 
		for (j = 0; j < N; j++) {
			comm_real[i + j * N] = HBAR_INV * (hamiltonian[j] * rho_imag[i + j * N] - hamiltonian[i] * rho_imag[i + j * N]);
			comm_imag[i + j * N] = HBAR_INV * (hamiltonian[i] * rho_real[i + j * N] - hamiltonian[j] * rho_real[i + j * N]);
		}
	}
}




void lindblad_operator(double * rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int SIZE) {
	int unsigned i, j, k, m, M, N;
	double rate, sum;	

	double* V;
	V = (double *) malloc(sizeof(double) * SIZE * SIZE);

	double *first_real, *second_real;//, *first_imag, *second_imag, *third_imag;
	double* helper;
	first_real  = (double *) malloc(sizeof(double) * SIZE * SIZE);
	second_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper      = (double *) malloc(sizeof(double) * SIZE * SIZE);



	#pragma acc data copyin(rho_imag[0:SIZE*SIZE], rho_real[0:SIZE*SIZE]) present_or_copyin(all_Vs[0:3*SIZE*SIZE*SIZE]) create(V[0:SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE])
	for (m = 1; m < SIZE - 1; m++) {

		#pragma acc loop independent
		for (M = 0; M < SIZE; M++) {
			#pragma acc loop independent
			for (N = 0; N < SIZE; N++) {

				#pragma acc kernels present(V[0:SIZE*SIZE])
				#pragma acc loop independent
				for (i = 0; i < SIZE; i++) {
					#pragma acc loop independent
					for (j = 0; j < SIZE; j++) {
						V[j + i * SIZE] = 0;
					}
				}
				V[M + N * SIZE] = all_Vs[M + N * SIZE + m * SIZE * SIZE + 1 * SIZE*SIZE*SIZE];
	

				#pragma acc kernels present(rho_real[0:SIZE*SIZE], V[0:SIZE*SIZE], helper[0:SIZE*SIZE])
				#pragma acc loop independent
				for (j = 0; j < SIZE; j++) {

					#pragma acc loop independent 
					for (i = 0; i < SIZE; i++) {
						helper[i + j * SIZE] = 0.;
						for (k = 0; k < SIZE; k++) {
							helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
						}
					}
				}
				//matrix_mul_real(helper, V_dagg, first_real, SIZE);
				// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
				#pragma acc kernels present(first_real[0:SIZE*SIZE], helper[0:SIZE*SIZE], V[0:SIZE*SIZE])
				#pragma acc loop independent 
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent 
					for (i = 0; i < SIZE; i++) {
						first_real[i + j * SIZE] = 0.;
						for (k = 0; k < SIZE; k++) {
							first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
						}
					}
				}
			//--------

				// getting the second part of the equation, V_dagg * V * rho / 2.
			//--------
				//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
				#pragma acc kernels present(second_real[0:SIZE*SIZE], V[0:SIZE*SIZE], helper[0:SIZE*SIZE])
				#pragma acc loop independent
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent 
					for (i = 0; i < SIZE; i++) {
						second_real[i + j * SIZE] = 0;
						for (k = 0; k < SIZE; k ++) {
							second_real[i + j * SIZE] += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
						}
					}
				}
			//--------

				// subtract first_real - second_real
				#pragma acc kernels present(first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE])
				#pragma acc loop independent 
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent 
					for (i = 0; i < SIZE; i++) {
						first_real[i + j * SIZE] = first_real[i + j * SIZE] - second_real[i + j * SIZE];
					}
				}


			//--------
				//matrix_mul_real(V_dagg, V, helper, SIZE);
				#pragma acc kernels present(helper[0:SIZE*SIZE], V[0:SIZE*SIZE])
				#pragma acc loop independent 
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent
					for (i = 0; i < SIZE; i++) {
						helper[i + j * SIZE] = 0;
						for (k = 0; k < SIZE; k ++) {
							helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
						}
					}
				}
				//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
				#pragma acc kernels present(helper[0:SIZE*SIZE], rho_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE]) 	
				#pragma acc loop independent 
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent 
					for (i = 0; i < SIZE; i++) {
						second_real[i + j * SIZE] = 0.;
						for (k = 0; k < SIZE; k ++) {
							second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
						}
					}
				}


			//--------
	

				// subtract first_real - second_real and multiply by the rate
				#pragma acc loop independent
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						lindblad_real[i + j * SIZE] += gammas[M + N * SIZE + m * SIZE * SIZE] * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
					}
				}

			}
		}
	}

	#pragma acc kernels present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3], eigVects[0:SIZE*SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], V[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copy(lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE]) 
	#pragma acc loop private(V[0:3*SIZE*SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE])
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 0 * SIZE*SIZE*SIZE];
			}
		}

		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0;
				for (k = 0; k < SIZE; k++) {
					helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
			}
		} 

		#pragma acc loop independent// worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop  independent// vector
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = 0;
				for (k = 0; k < SIZE; k++) {
					first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
				}
			}
		}

	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0;
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}

	//--------

		// subtract first_real - second_real
		#pragma acc loop independent
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = first_real[i + j * SIZE] - second_real[i + j * SIZE];
			}
		}


		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
		#pragma acc loop independent// worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop  independent //vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k ++) {
					helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}
	//--------

		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i + j * SIZE] += rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}		

	}

	// decay of excitons into the target (reaction center) state
        
	#pragma acc kernels present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3], eigVects[0:SIZE*SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], V[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copy(lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
    #pragma acc loop private(V[0:3*SIZE*SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE])
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];
		#pragma acc loop independent //worker 
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent //vector
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 2 * SIZE*SIZE*SIZE];
			}
		}

		// getting the first part of the equation, V * rho * V_dagg
	//--------
		//matrix_mul_real(V, rho_real, helper, SIZE);
		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector 
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k++) {
					helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(helper, V_dagg, first_real, SIZE);
		// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
		#pragma acc loop independent// worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = 0;
				for (k = 0; k < SIZE; k++) {
					first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
				}
			}
		}
	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}
		//--------

		// subtract first_real - second_real
		#pragma acc loop independent //worker 
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = first_real[i + j * SIZE] - second_real[i + j * SIZE];
			}
		}


		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
		#pragma acc loop independent// worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent// vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k ++) {
					helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent //worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent //vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}
	//--------

		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent //worker 
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i + j * SIZE] += rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}
	} // end of acc data directive

	free((void*) V);
	free((void*) first_real);
	free((void*) second_real);
	free((void*) helper);


}



void get_density_update(double *rho_real, double *rho_imag, double *energies, double *comm_real, double *comm_imag, 
	                    double *gammas, double *eigvects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int N) {

	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

	int unsigned i, j;
	#pragma acc kernels copy(lindblad_real[0:N*N], lindblad_imag[0:N*N])
	#pragma acc loop parallel independent
	for (i = 0; i < N; i++) {
		#pragma acc loop parallel independent
		for (j = 0; j < N; j++) {
			lindblad_real[i + j * N] = 0.;
			lindblad_imag[i + j * N] = 0.;
		}
	}

	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, N);

}

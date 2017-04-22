
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

void hamiltonian_commutator(double **rho_real, double **rho_imag, double *hamiltonian, double **comm_real, double **comm_imag, int N) {
	// as dense as possible implementation of the hamiltonian commutator -i * [H, rho] / hbar
	int unsigned i, j;

	#pragma acc kernels present(hamiltonian[0:N])
	#pragma acc loop independent collapse(2)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			comm_real[i][j] = (hamiltonian[j] * rho_imag[i][j] - hamiltonian[i] * rho_imag[i][j]) * HBAR_INV;
			comm_imag[i][j] = (hamiltonian[i] * rho_real[i][j] - hamiltonian[j] * rho_real[i][j]) * HBAR_INV;

		}
	}

}




void lindblad_operator(double **rho_real, double **rho_imag, double ***gammas, double **eigVects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, int SIZE) {
	int unsigned i, j, k, m, M, N;
	double rate, sum;	

	double *V[SIZE];
	for (i = 0; i < SIZE; i++) {
		V[i] = (double *) malloc(sizeof(double) * SIZE);
	}

	double *first_real[SIZE], *second_real[SIZE], *helper[SIZE];
	for (i = 0; i < SIZE; i++) {
		first_real[i] = (double *) malloc(sizeof(double) * SIZE);
		second_real[i] = (double *) malloc(sizeof(double) * SIZE);
		helper[i] = (double *) malloc(sizeof(double) * SIZE);
	}


//	printf("starting liouville\n");

 
//	#pragma acc kernels present(gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][3])    create(V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE])
	#pragma acc data copyin(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])    present(gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][3])
	for (m = 1; m < SIZE - 1; m++) {

		#pragma acc loop independent 
		for (M = 0; M < SIZE; M++) {
			#pragma acc loop independent
			for (N = 0; N < SIZE; N++) {

				#pragma acc kernels create(V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE])
				{

				#pragma acc loop independent collapse(2)
				for (i = 0; i < SIZE; i++) {
					for (j = 0; j < SIZE; j++) {
						V[j][i] = 0;
					}
				}
				V[M][N] = all_Vs[M][N][m][1];


				#pragma acc loop independent worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent vector
					for (i = 0; i < SIZE; i++) {
						helper[i][j] = 0.;
						#pragma acc loop independent
						for (k = 0; k < SIZE; k++) {
							helper[i][j] += V[i][k] * rho_real[k][j];
						}
					}
				}
	

				//matrix_mul_real(helper, V_dagg, first_real, SIZE);
				// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
				#pragma acc loop independent worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent vector
					for (i = 0; i < SIZE; i++) {
						first_real[i][j] = 0.;
						#pragma acc loop independent
						for (k = 0; k < SIZE; k++) {
							first_real[i][j] += helper[i][k] * V[j][k];
						}
					}
				}
			//--------

				// getting the second part of the equation, V_dagg * V * rho / 2.
			//--------
				//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
				#pragma acc loop independent worker
				for (j = 0; j < SIZE; j++) { 
					#pragma acc loop independent vector
					for (i = 0; i < SIZE; i++) {
						second_real[i][j] = 0;
						#pragma acc loop independent
						for (k = 0; k < SIZE; k ++) {
							second_real[i][j] += V[k][i] * helper[k][j] / 2.;
						}
					}
				}
			//--------

				// subtract first_real - second_real
				#pragma acc loop independent collapse(2)
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						first_real[i][j] = first_real[i][j] - second_real[i][j];
					}
				}


			//--------
				//matrix_mul_real(V_dagg, V, helper, SIZE);

				#pragma acc loop independent worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent vector
					for (i = 0; i < SIZE; i++) {
						helper[i][j] = 0;
						#pragma acc loop independent 
						for (k = 0; k < SIZE; k ++) {
							helper[i][j] += V[k][i] * V[k][j];
						}
					}
				}
				//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
				#pragma acc loop independent worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent vector 
					for (i = 0; i < SIZE; i++) {
						second_real[i][j] = 0.;
						#pragma acc loop independent
						for (k = 0; k < SIZE; k ++) {
							second_real[i][j] += rho_real[i][k] * helper[k][j] / 2.;
						}
					}
				}


			//--------
	

				// subtract first_real - second_real and multiply by the rate
				#pragma acc loop independent collapse(2)
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						lindblad_real[i][j] += gammas[M][N][m] * (first_real[i][j] - second_real[i][j]);
					}
				}
			
				} // end kernels			

			}
		}
	}


//	#pragma acc kernels present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3], eigVects[0:SIZE*SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], V[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copy(lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE]) 
//	#pragma acc loop private(V[0:3*SIZE*SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE])
        #pragma acc kernels present(links_to_loss[0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][3])    create(V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE])
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				V[j][i] = all_Vs[j][i][m][0];
			}
		}


		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				helper[i][j] = 0;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k++) {
					helper[i][j] += V[i][k] * rho_real[k][j];
				}
			}
		} 

		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop  independent vector
			for (i = 0; i < SIZE; i++) {
				first_real[i][j] = 0;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k++) {
					first_real[i][j] += helper[i][k] * V[j][k];
				}
			}
		}

	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				second_real[i][j] = 0;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					second_real[i][j] += V[k][i] * helper[k][j] / 2.;
				}
			}
		}

	//--------

		// subtract first_real - second_real
		#pragma acc loop independent collapse(2)
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				first_real[i][j] = first_real[i][j] - second_real[i][j];
			}
		}

		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop  independent vector
			for (i = 0; i < SIZE; i++) {
				helper[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					helper[i][j] += V[k][i] * V[k][j];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				second_real[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					second_real[i][j] += rho_real[i][k] * helper[k][j] / 2.;
				}
			}
		}
	//--------

		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent collapse(2)
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i][j] += rate * (first_real[i][j] - second_real[i][j]);
			}
		}		

	}


	// decay of excitons into the target (reaction center) state
        
//	#pragma acc kernels present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3], eigVects[0:SIZE*SIZE], links_to_loss[0:SIZE], links_to_target[0:SIZE], V[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE]) copy(lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
//    #pragma acc loop private(V[0:3*SIZE*SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], helper[0:SIZE*SIZE])


	#pragma acc kernels present(links_to_target[0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][3])     create(V[0:SIZE][0:SIZE], first_real[0:SIZE][0:SIZE], second_real[0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE])
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];

		#pragma acc loop independent collapse(2)
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				V[j][i] = all_Vs[j][i][m][2];
			}
		}


		// getting the first part of the equation, V * rho * V_dagg
	//--------
		//matrix_mul_real(V, rho_real, helper, SIZE);
		#pragma acc loop independent worker 
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				helper[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k++) {
					helper[i][j] += V[i][k] * rho_real[k][j];
				}
			}
		}

		//matrix_mul_real(helper, V_dagg, first_real, SIZE);
		// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				first_real[i][j] = 0;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k++) {
					first_real[i][j] += helper[i][k] * V[j][k];
				}
			}
		}
	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				second_real[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					second_real[i][j] += V[k][i] * helper[k][j] / 2.;
				}
			}
		}
		//--------

		// subtract first_real - second_real
		#pragma acc loop independent collapse(2)
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				first_real[i][j] = first_real[i][j] - second_real[i][j];
			}
		}


		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				helper[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					helper[i][j] += V[k][i] * V[k][j];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop independent worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent vector
			for (i = 0; i < SIZE; i++) {
				second_real[i][j] = 0.;
				#pragma acc loop independent
				for (k = 0; k < SIZE; k ++) {
					second_real[i][j] += rho_real[i][k] * helper[k][j] / 2.;
				}
			}
		}
	//--------

		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent collapse(2)
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i][j] += rate * (first_real[i][j] - second_real[i][j]);
			}
		}

	} // end of acc data directive


}


void get_density_update(double **rho_real, double **rho_imag, double *energies, double **comm_real, double **comm_imag, 
	                    double ***gammas, double **eigvects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, int N) {

	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

	int unsigned i, j;
	#pragma acc kernels
	#pragma acc loop parallel independent collapse(2)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			lindblad_real[i][j] = 0.;
			lindblad_imag[i][j] = 0.;
		}
	}

	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, N);

}

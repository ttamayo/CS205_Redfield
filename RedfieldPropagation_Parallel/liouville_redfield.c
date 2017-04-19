
#include <stdio.h>
#include <stdlib.h>

#include "headers.h"
<<<<<<< HEAD
//#include "multiple_matrix_operations.h"
=======
<<<<<<< HEAD
//#include "multiple_matrix_operations.h"
=======
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

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

<<<<<<< HEAD

/********************************************************************/

#pragma acc routine gang
void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N) {
=======
/********************************************************************/

<<<<<<< HEAD
#pragma acc routine gang
void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N) {
	// as dense as possible implementation of the hamiltonian commutator -i * [H, rho] / hbar
	int unsigned i, j;
	#pragma acc data present_or_copyin(hamiltonian[0:N]) copyin(rho_real[0:N*N], rho_imag[0:N*N]) present_or_copy(comm_real[0:N*N], comm_imag[0:N*N])
	#pragma acc loop independent
=======
 void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N) {
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
 	// FIXME
 	// the hamiltonian is a real valued diagonal matrix
 	// ... but we don't care for now and use simple matrix multiplication

<<<<<<< HEAD
	int unsigned i, j;
	// hamiltonian is a simple one dimensional array
	
	double *help1_real, *help1_imag, *help2_real, *help2_imag;
	help1_real = (double *) malloc(sizeof(double) * (N * N));
=======
 	double *h_real, *h_imag; 
 	h_real = (double *) malloc(sizeof(double) * (N * N));
 	h_imag = (double *) malloc(sizeof(double) * (N * N));
 	gen_zero_matrix_real(h_real, N);
 	gen_zero_matrix_real(h_imag, N);

 	double *help1_real, *help1_imag, *help2_real, *help2_imag;
 	help1_real = (double *) malloc(sizeof(double) * (N * N));
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
	help2_real = (double *) malloc(sizeof(double) * (N * N));
 	help1_imag = (double *) malloc(sizeof(double) * (N * N));
	help2_imag = (double *) malloc(sizeof(double) * (N * N));


	// get the first part of the commutator
<<<<<<< HEAD
	#pragma acc data present(rho_real[0:N*N], rho_imag[0:N*N], hamiltonian[0:N]) create(help1_real[0:N*N], help1_imag[0:N*N]) 
	#pragma acc loop gang
	for (i = 0; i < N; i++) {
		#pragma acc loop gang, worker
		for (j = 0; j < N; j++) {
			help1_real[i + j * N] = hamiltonian[j] * rho_real[i + j * N];
			help1_imag[i + j * N] = 0;
		}
	}

	// get the second part of the commutator rho * H
	#pragma acc data present(rho_real[0:N*N], rho_imag[0:N*N], hamiltonian[0:N]) create(help2_real[0:N*N], help2_imag[0:N*N])
	#pragma acc loop gang
	for (i = 0; i < N; i++) {
		#pragma acc loop gang, worker
		for (j = 0; j < N; j++) {
			help2_real[i + j * N] = hamiltonian[i] * rho_real[i + j * N];
			help2_imag[i + j * N] = 0;
		}
	}
	
	
	// combine the two parts
	#pragma acc data present(help1_real[0:N*N], help1_imag[0:N*N], help2_real[0:N*N], help2_imag[0:N*N], comm_real[0:N*N], comm_imag[0:N*N])
	matrix_sub_real(help1_imag, help2_imag, comm_real, N);
	matrix_sub_real(help2_real, help1_real, comm_imag, N);

	matrix_mul_scalar(comm_real, HBAR_INV, N);
	matrix_mul_scalar(comm_imag, HBAR_INV, N);
	#pragma end data

=======
	int unsigned i;
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	for (i = 0; i < N; i++) {
		#pragma acc loop independent
		for (j = 0; j < N; j++) {
			comm_real[i + j * N] = HBAR_INV * (hamiltonian[j] * rho_imag[i + j * N] - hamiltonian[i] * rho_imag[i + j * N]);
			comm_imag[i + j * N] = HBAR_INV * (hamiltonian[i] * rho_real[i + j * N] - hamiltonian[j] * rho_real[i + j * N]);
		}
	}
}


<<<<<<< HEAD
=======
	free((void*) h_real);
	free((void*) h_imag);
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
	free((void*) help1_real);
	free((void*) help1_imag);
	free((void*) help2_real);
	free((void*) help2_imag);
	
<<<<<<< HEAD
=======
	#pragma end data
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
	
 }
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa



#pragma acc routine gang
<<<<<<< HEAD
void lindblad_operator(double * restrict rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int SIZE) {
	int unsigned i, j, k, m, M, N;
	double rate, sum;	
=======
void lindblad_operator(double *rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, int SIZE) {
	int unsigned m, M, N;
<<<<<<< HEAD
	double restrict rate;
=======
	double rate;
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

	double * restrict V;
	V = (double *) malloc(sizeof(double) * SIZE * SIZE);

	double *first_real, *second_real;//, *first_imag, *second_imag, *third_imag;
	double * restrict helper;
	first_real  = (double *) malloc(sizeof(double) * SIZE * SIZE);
	second_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper      = (double *) malloc(sizeof(double) * SIZE * SIZE);


<<<<<<< HEAD
	{
//	#pragma acc loop independent gang

	#pragma acc data present_or_copyin(gammas[0:SIZE*SIZE*SIZE], all_Vs[0:SIZE*SIZE*SIZE*3]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE])
	#pragma acc loop parallel gang
	for (m = 1; m < SIZE - 1; m++) {
		// get the V matrix for this particular transition		
		#pragma acc loop parallel independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop parallel independent
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 1 * SIZE*SIZE*SIZE];
			}
		}

		#pragma acc loop parallel worker independent
		for (M = 0; M < SIZE; M++) {
			#pragma acc loop parallel vector independent
			for (N = 0; N < SIZE; N++) {

				rate = gammas[M + N * SIZE + m * SIZE * SIZE];
		
				// getting the first part of the equation, V * rho * V_dagg
			//--------
				//matrix_mul_real(V, rho_real, helper, SIZE);
				#pragma acc loop worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector
					for (i = 0; i < SIZE; i++) {
						helper[i + j * SIZE] = 0.;
						#pragma acc loop independent seq
						for (k = 0; k < SIZE; k++) {
							helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
						}
					}
				}
				//matrix_mul_real(helper, V_dagg, first_real, SIZE);
				// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
				#pragma acc loop worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector
					for (i = 0; i < SIZE; i++) {
						first_real[i + j * SIZE] = 0.;
						#pragma acc loop seq
						for (k = 0; k < SIZE; k++) {
							first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
						}
					}
				}
			//--------

				// getting the second part of the equation, V_dagg * V * rho / 2.
			//--------
				//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
				#pragma acc loop worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector
					for (i = 0; i < SIZE; i++) {
						second_real[i + j * SIZE] = 0;
						#pragma acc loop seq
						for (k = 0; k < SIZE; k ++) {
							second_real[i + j * SIZE] += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
						}
					}
				}
			//--------

				// subtract first_real - second_real
				#pragma acc loop worker independent
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector independent
					for (i = 0; i < SIZE; i++) {
						first_real[i + j * SIZE] = first_real[i + j * SIZE] - second_real[i + j * SIZE];
					}
				}


			//--------
				//matrix_mul_real(V_dagg, V, helper, SIZE);
				#pragma acc loop worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector
					for (i = 0; i < SIZE; i++) {
						helper[i + j * SIZE] = 0;
						#pragma acc loop seq
						for (k = 0; k < SIZE; k ++) {
							helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
						}
					}
				}
				//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
				#pragma acc loop worker
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop vector
					for (i = 0; i < SIZE; i++) {
						second_real[i + j * SIZE] = 0.;
						#pragma acc loop seq	
						for (k = 0; k < SIZE; k ++) {
							second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
						}
					}
				}
			//--------

				// subtract first_real - second_real and multiply by the rate
				#pragma acc loop independent
				for (j = 0; j < SIZE; j++) {
					#pragma acc loop independent
					for (i = 0; i < SIZE; i++) {
						lindblad_real[i + j * SIZE] += rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
					}
				}
=======
//       #pragma acc data copyin(rho_real[0:N*N], rho_imag[0:N*N], h_real[0:N*N], h_imag[0:N*N]) create(help1_real[0:N*N], help1_imag[0:N*N], help2_real[0:N*N], help2_imag[0:N*N]) copy(comm_real[0:N*N], comm_imag[0:N*N])
 
	// decay of excitons into other excitonic states
//	#pragma acc data copyin(gammas[0:SIZE*SIZE*SIZE])
//	#pragma acc parallel loop tile(16, 16) gang vector
<<<<<<< HEAD
//	#pragma acc data copyin(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) copy(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE])
	#pragma acc declare copyin(SIZE) 
	#pragma acc data create(helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) present_or_copyin(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE]) present(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
//	#pragma acc loop gang
=======
	#pragma acc data copyin(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) copy(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE])
 

>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
	for (m = 1; m < SIZE - 1; m++) {
		get_V(V, eigVects, m, m, SIZE);
		get_V(V_dagg, eigVects, m, m, SIZE);
		transpose(V_dagg, SIZE);
<<<<<<< HEAD

//		#pragma acc loop gang, worker
		for (M = 0; M < SIZE; M++) {
//			#pragma acc loop gang, worker
=======
		#pragma acc loop independent
		for (M = 0; M < SIZE; M++) {
			#pragma acc loop independent
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
			for (N = 0; N < SIZE; N++) {

				rate = gammas[M + N * SIZE + m * SIZE * SIZE];
				rate = (double)(rate);
<<<<<<< HEAD
		
=======

>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
//				get_V(V, eigVects, m, m, SIZE);
//				get_V(V_dagg, eigVects, m, m, SIZE);
//				transpose(V_dagg, SIZE);

				matrix_mul_real(V, rho_real, helper, SIZE);
				matrix_mul_real(helper, V_dagg, first_real, SIZE);

				matrix_mul_real(V, rho_real, helper, SIZE);
				matrix_mul_real(V_dagg, helper, second_real, SIZE);
				matrix_mul_scalar(second_real, 0.5, SIZE);

				matrix_mul_real(V_dagg, V, helper, SIZE);
				matrix_mul_real(rho_real, helper, third_real, SIZE);
				matrix_mul_scalar(third_real, 0.5, SIZE);	

				matrix_sub_real(first_real, second_real, first_real, SIZE);
				matrix_sub_real(first_real, third_real, first_real, SIZE);
				matrix_mul_scalar(first_real, rate, SIZE);

>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa


			}
		}
	}

<<<<<<< HEAD

	#pragma acc data present_or_copyin(all_Vs[0:SIZE*SIZE*SIZE*3], V[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE])
	#pragma acc loop gang
=======
<<<<<<< HEAD
//	#pragma end data

	// decay of excitons into the ground (loss) state
	#pragma acc data create(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) present(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
	#pragma acc loop gang
=======
	#pragma end data

	// decay of excitons into the ground (loss) state
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 0 * SIZE*SIZE*SIZE];
			}
		}

		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k++) {
					helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
			}
		} 

		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = 0;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k++) {
					first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
				}
			}
		}
	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0;
				#pragma acc loop seq
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
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k ++) {
					helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}
	//--------

		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i + j * SIZE] += rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}		
	}

<<<<<<< HEAD

	// decay of excitons into the target (reaction center) state
        
	#pragma acc data present_or_copyin(all_Vs[0:SIZE*SIZE*SIZE*3], V[0:SIZE*SIZE], helper[0:SIZE*SIZE]) copyin(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE])
	#pragma acc loop gang
=======
<<<<<<< HEAD
///	#pragma end data

	// decay of excitons into the target (reaction center) state
//	#pragma acc data create(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) present(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
	#pragma acc data create(V[0:SIZE*SIZE], V_dagg[0:SIZE*SIZE], helper[0:SIZE*SIZE], first_real[0:SIZE*SIZE], second_real[0:SIZE*SIZE], third_real[0:SIZE*SIZE]) present(rho_real[0:SIZE*SIZE], rho_imag[0:SIZE*SIZE], gammas[0:SIZE*SIZE*SIZE], lindblad_real[0:SIZE*SIZE], lindblad_imag[0:SIZE*SIZE])
	#pragma acc loop gang
=======
	// decay of excitons into the target (reaction center) state
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];

		#pragma acc loop independent
		for (i = 0; i < SIZE; i++) {
			#pragma acc loop independent
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 2 * SIZE*SIZE*SIZE];
			}
		}

		// getting the first part of the equation, V * rho * V_dagg
	//--------
		//matrix_mul_real(V, rho_real, helper, SIZE);
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop gang vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k++) {
					helper[i + j * SIZE] += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(helper, V_dagg, first_real, SIZE);
		// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				first_real[i + j * SIZE] = 0;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k++) {
					first_real[i + j * SIZE] += helper[i + k * SIZE] * V[j + k * SIZE];
				}
			}
		}
	//--------

		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				#pragma acc loop seq
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
		#pragma acc loop worker
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				helper[i + j * SIZE] = 0.;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k ++) {
					helper[i + j * SIZE] += V[k + i * SIZE] * V[k + j * SIZE];
				}
			}
		}

		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
		#pragma acc loop  
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop vector
			for (i = 0; i < SIZE; i++) {
				second_real[i + j * SIZE] = 0.;
				#pragma acc loop seq
				for (k = 0; k < SIZE; k ++) {
					second_real[i + j * SIZE] += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
			}
		}
	//--------

<<<<<<< HEAD
		// subtract first_real - second_real and multiply by the rate
		#pragma acc loop independent
		for (j = 0; j < SIZE; j++) {
			#pragma acc loop independent
			for (i = 0; i < SIZE; i++) {
				lindblad_real[i + j * SIZE] += rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}
	} // end of acc data directive
	} 
=======
<<<<<<< HEAD
	#pragma end data 
=======
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

	free((void*) V);
	free((void*) first_real);
	free((void*) second_real);
	free((void*) helper);


}



<<<<<<< HEAD
#pragma acc routine gang
=======
<<<<<<< HEAD
#pragma acc routine gang
=======
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
void get_density_update(double *rho_real, double *rho_imag, double *energies, double *comm_real, double *comm_imag, 
	                    double *gammas, double *eigvects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int N) {

	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

<<<<<<< HEAD
	int unsigned i, j;
	#pragma acc loop independent
	for (i = 0; i < N; i++) {
		#pragma acc loop independent
		for (j = 0; j < N; j++) {
			lindblad_real[i + j * N] = 0.;
			lindblad_imag[i + j * N] = 0.;
		}
	}

	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, N);

}
=======
<<<<<<< HEAD
}
=======
}
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

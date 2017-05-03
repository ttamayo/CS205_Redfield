
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

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

#define BLOCK_SIZE 16

/********************************************************************/

void hamiltonian_commutator(double * rho_real, double * rho_imag, double * hamiltonian, double * comm_real, double * comm_imag, int N) {
	// as dense as possible implementation of the hamiltonian commutator -i * [H, rho] / hbar
	int unsigned i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			comm_real[i + j * N] = (hamiltonian[j] * rho_imag[i + j * N] - hamiltonian[i] * rho_imag[i + j * N]) * HBAR_INV;
			comm_imag[i + j * N] = (hamiltonian[i] * rho_real[i + j * N] - hamiltonian[j] * rho_real[i + j * N]) * HBAR_INV;
		}
	}
}




void lindblad_operator(double * rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int SIZE) {
	int unsigned i, j, k, m, M, N;
	double rate, temp;	

	double* V;
	V = (double *) malloc(sizeof(double) * SIZE * SIZE);

	double *first_real, *second_real;//, *first_imag, *second_imag, *third_imag;
	double* helper;
	first_real  = (double *) malloc(sizeof(double) * SIZE * SIZE);
	second_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper      = (double *) malloc(sizeof(double) * SIZE * SIZE);

	
	// for tiling
	int unsigned ii,jj,index,jndex;
	int bsize = BLOCK_SIZE;


	for (m = 1; m < SIZE - 1; m++) {

		for (M = 0; M < SIZE; M++) {

			for (N = 0; N < SIZE; N++) {

				for (i = 0; i < SIZE; i++) {
					for (j = 0; j < SIZE; j++) {
						V[j + i * SIZE] = 0;
					}
				}

				V[M + N * SIZE] = all_Vs[M + N * SIZE + m * SIZE * SIZE + 1 * SIZE*SIZE*SIZE];

				// (1)
/*
				#pragma omp parallel shared(helper,V,rho_real,SIZE,temp) private(j,i,k) 
				#pragma omp for nowait
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = 0.;
						for (k = 0; k < SIZE; k++) {
							temp += V[i + k * SIZE] * rho_real[k + j * SIZE];
						}
						helper[i+j*SIZE] = temp;
					}
				}
*/
				#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
				for (ii = 0; ii < SIZE; ii += bsize) {
					for (jj = 0; jj < SIZE; jj += bsize) {
						for (i = 0; i < SIZE; i++ ) {
							for (j = 0; j < SIZE; j++) {
								index = ii + i;
								jndex = jj + j;
								if (index < SIZE && jndex < SIZE) {  
									temp = 0;
									for (k = 0; k < SIZE; k++ ) {
										temp += V[index+k*SIZE] * rho_real[k+jndex*SIZE];
									}
									helper[i+j*SIZE] = temp;
								}
							}
						}
					}
				}





				// (2)
				//matrix_mul_real(helper, V_dagg, first_real, SIZE);
				// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
/*	
				#pragma omp parallel shared(first_real,helper,V,SIZE,temp) private(j,i,k) 
				#pragma omp for nowait
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = 0.;
						for (k = 0; k < SIZE; k++) {
							temp += helper[i + k * SIZE] * V[j + k * SIZE];
						}
						first_real[i+j*SIZE] = temp;
					}
				}
*/	
				#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                                for (ii = 0; ii < SIZE; ii += bsize) {
                                        for (jj = 0; jj < SIZE; jj += bsize) {
                                                for (i = 0; i < SIZE; i++ ) {
                                                        for (j = 0; j < SIZE; j++) {
                                                                index = ii + i;
                                                                jndex = jj + j;
                                                                if (index < SIZE && jndex < SIZE) {
                                                                        temp = 0;
                                                                        for (k = 0; k < SIZE; k++ ) {
                                                                                temp += helper[index+k*SIZE] * V[jndex+k*SIZE];
                                                                        }
                                                                        first_real[i+j*SIZE] = temp;
                                                                }
                                                        }
                                                }
                                        }
                                }



			// (3)
			//--------

				// getting the second part of the equation, V_dagg * V * rho / 2.
			//--------
				//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
/*				
				#pragma omp parallel shared(second_real,V,helper,SIZE,temp) private(j,i,k)
				#pragma omp for nowait
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = 0.;
						for (k = 0; k < SIZE; k ++) {
							temp += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
						}
						second_real[i+j*SIZE] = temp;
					}
				}
*/
				#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                                for (ii = 0; ii < SIZE; ii += bsize) {
                                        for (jj = 0; jj < SIZE; jj += bsize) {
                                                for (i = 0; i < SIZE; i++ ) {
                                                        for (j = 0; j < SIZE; j++) {
                                                                index = ii + i;
                                                                jndex = jj + j;
                                                                if (index < SIZE && jndex < SIZE) {
                                                                        temp = 0;
                                                                        for (k = 0; k < SIZE; k++ ) {
                                                                                temp += V[k+index*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                                        }
                                                                        second_real[i+j*SIZE] = temp;
                                                                }
                                                        }
                                                }
                                        }
                                }
				

			//--------

				// subtract first_real - second_real
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = first_real[i+j*SIZE];
						first_real[i + j * SIZE] = temp - second_real[i + j * SIZE];
					}
				}


			//-------
			// (4)
				//matrix_mul_real(V_dagg, V, helper, SIZE);
/*				
				#pragma omp parallel shared(helper,V,SIZE,temp) private(j,i,k) 
				#pragma omp for nowait
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = 0.;
						for (k = 0; k < SIZE; k ++) {
							temp += V[k + i * SIZE] * V[k + j * SIZE];
						}
						helper[i+j*SIZE] = temp;
					}
				}
*/
				#pragma omp parallel for  private(ii,jj,i,j,k) schedule(auto) collapse(2)
                                for (ii = 0; ii < SIZE; ii += bsize) {
                                        for (jj = 0; jj < SIZE; jj += bsize) {
                                                for (i = 0; i < SIZE; i++ ) {
                                                        for (j = 0; j < SIZE; j++) {
                                                                index = ii + i;
                                                                jndex = jj + j;
                                                                if (index < SIZE && jndex < SIZE) {
                                                                        temp = 0;
                                                                        for (k = 0; k < SIZE; k++ ) {
                                                                                temp += V[k+index*SIZE] * V[k+jndex*SIZE];
                                                                        }
                                                                        helper[i+j*SIZE] = temp;
                                                                }
                                                        }
                                                }
                                        }
                                }


				// (5)
				//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
				
/*
				#pragma omp parallel shared(second_real,rho_real,helper,SIZE,temp) private(j,i,k)
				#pragma omp for nowait
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						temp = 0.;
						for (k = 0; k < SIZE; k ++) {
							temp += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
						}
						second_real[i+j*SIZE] = temp;
					}
				}
*/
				#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                                for (ii = 0; ii < SIZE; ii += bsize) {
                                        for (jj = 0; jj < SIZE; jj += bsize) {
                                                for (i = 0; i < SIZE; i++ ) {
                                                        for (j = 0; j < SIZE; j++) {
                                                                index = ii + i;
                                                                jndex = jj + j;
                                                                if (index < SIZE && jndex < SIZE) {
                                                                        temp = 0;
                                                                        for (k = 0; k < SIZE; k++ ) {
                                                                                temp += rho_real[index+k*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                                        }
                                                                        second_real[i+j*SIZE] = temp;
                                                                }
                                                        }
                                                }
                                        }
                                }



			//--------
	

				// subtract first_real - second_real and multiply by the rate
				for (j = 0; j < SIZE; j++) {
					for (i = 0; i < SIZE; i++) {
						lindblad_real[i + j * SIZE] += gammas[M + N * SIZE + m * SIZE * SIZE] * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
					}
				}

			}
		}
	}


	//print_matrix_real(lindblad_real, SIZE);



	for (m = 0; m < SIZE; m++) {
		

		rate = links_to_loss[m];

		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 0 * SIZE*SIZE*SIZE];
			}
		}

		// (6)
/*
		#pragma omp parallel shared(helper,V,rho_real,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k++) {
					temp += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
				helper[i+j*SIZE] = temp;
			}
		} 
*/

		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                	for (jj = 0; jj < SIZE; jj += bsize) {
                        	for (i = 0; i < SIZE; i++ ) {
                                	for (j = 0; j < SIZE; j++) {
                                        	index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                	temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                        	temp += V[index+k*SIZE] * rho_real[k+jndex*SIZE];
                                                        }
                                                        helper[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }	

		// (7)
/*
		#pragma omp parallel shared(first_real,helper,V,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0;
				for (k = 0; k < SIZE; k++) {
					temp += helper[i + k * SIZE] * V[j + k * SIZE];
				}
				first_real[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += helper[index+k*SIZE] * V[jndex+k*SIZE];
                                                        }
                                                        first_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }





	//--------
		// (8)
		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
/*
		#pragma omp parallel shared(second_real,V,helper,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0;
				for (k = 0; k < SIZE; k ++) {
					temp += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
				}
				second_real[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += V[k+index*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                        }
                                                        second_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }





	//--------
		// subtract first_real - second_real
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = first_real[i+j*SIZE];
				first_real[i + j * SIZE] = temp - second_real[i + j * SIZE];
			}
		}

		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
		// (9)
/*
		#pragma omp parallel shared(helper,V,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k ++) {
					temp += V[k + i * SIZE] * V[k + j * SIZE];
				}
				helper[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += V[k+index*SIZE] * V[k+jndex*SIZE];
                                                        }
                                                        helper[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }


		// (10)
		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
/*
		#pragma omp parallel shared(second_real,rho_real,helper,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k ++) {
					temp += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
				second_real[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += rho_real[index+k*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                        }
                                                        second_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }



	//--------

		// subtract first_real - second_real and multiply by the rate
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = lindblad_real[i+j*SIZE];
				lindblad_real[i + j * SIZE] = temp + rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}		

	}

	// decay of excitons into the target (reaction center) state
        
	for (m = 0; m < SIZE; m++) {
		
		rate = links_to_target[m];
		
		for (i = 0; i < SIZE; i++) {
			for (j = 0; j < SIZE; j++) {
				V[j + i * SIZE] = all_Vs[j + i * SIZE + m * SIZE*SIZE + 2 * SIZE*SIZE*SIZE];
			}
		}


		// getting the first part of the equation, V * rho * V_dagg
	//--------
	// (11)
		//matrix_mul_real(V, rho_real, helper, SIZE);
/*
		#pragma omp parallel shared(helper,V,rho_real,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k++) {
					temp += V[i + k * SIZE] * rho_real[k + j * SIZE];
				}
				helper[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += V[index+k*SIZE] * rho_real[k+jndex*SIZE];
                                                        }
                                                        helper[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }



		// (12)
		//matrix_mul_real(helper, V_dagg, first_real, SIZE);
		// V_dagg is 'computed implicitly' by using the correct indexing in the matrix multiplication
/*
		#pragma omp parallel shared(first_real,helper,V,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0;
				for (k = 0; k < SIZE; k++) {
					temp += helper[i + k * SIZE] * V[j + k * SIZE];
				}
				first_real[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += helper[index+k*SIZE] * V[jndex+k*SIZE];
                                                        }
                                                        first_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }



	//--------
		// (13)
		// getting the second part of the equation, V_dagg * V * rho / 2.
	//--------
		//matrix_mul_real(V_dagg, helper, second_real, SIZE) / 2.;
/*		
		#pragma omp parallel shared(second_real,V,helper,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k ++) {
					temp += V[k + i * SIZE] * helper[k + j * SIZE] / 2.;
				}
				second_real[i+j*SIZE] = temp;
			}
		}
*/	
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += V[k+index*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                        }
                                                        second_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }

		//--------

		// subtract first_real - second_real
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = first_real[i+j*SIZE];
				first_real[i + j * SIZE] = temp - second_real[i + j * SIZE];
			}
		}


		// (14)
		// getting the third part of the equation, rho * V_dagg * V / 2.
	//--------
		//matrix_mul_real(V_dagg, V, helper, SIZE);
/*
		#pragma omp parallel shared(helper,V,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k ++) {
					temp += V[k + i * SIZE] * V[k + j * SIZE];
				}
				helper[i+j*SIZE] = temp;
			}
		}
*/
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += V[k+index*SIZE] * V[k+jndex*SIZE];
                                                        }
                                                        helper[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }



		// (15)
		//matrix_mul_real(rho, helper, second_real, SIZE) / 2.;
/*
		#pragma omp parallel shared(second_real,rho_real,helper,SIZE,temp) private(j,i,k)
		#pragma omp for nowait
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = 0.;
				for (k = 0; k < SIZE; k ++) {
					temp += rho_real[i + k * SIZE] * helper[k + j * SIZE] / 2.;
				}
				second_real[i+j*SIZE] = temp;
			}
		}
*/		
		#pragma omp parallel for private(ii,jj,i,j,k) schedule(auto) collapse(2)
                for (ii = 0; ii < SIZE; ii += bsize) {
                        for (jj = 0; jj < SIZE; jj += bsize) {
                                for (i = 0; i < SIZE; i++ ) {
                                        for (j = 0; j < SIZE; j++) {
                                                index = ii + i;
                                                jndex = jj + j;
                                                if (index < SIZE && jndex < SIZE) {
                                                        temp = 0;
                                                        for (k = 0; k < SIZE; k++ ) {
                                                                temp += rho_real[index+k*SIZE] * helper[k+jndex*SIZE] * 0.5;
                                                        }
                                                        second_real[i+j*SIZE] = temp;
                                                }
                                        }
                                }
                        }
                }



	//--------

		// subtract first_real - second_real and multiply by the rate
		for (j = 0; j < SIZE; j++) {
			for (i = 0; i < SIZE; i++) {
				temp = lindblad_real[i+j*SIZE];
				lindblad_real[i + j * SIZE] = temp + rate * (first_real[i + j * SIZE] - second_real[i + j * SIZE]);
			}
		}

	} 

	free((void*) V);
	free((void*) first_real);
	free((void*) second_real);
	free((void*) helper);


}



void get_density_update(double *rho_real, double *rho_imag, double *energies, double *comm_real, double *comm_imag, 
	                    double *gammas, double *eigvects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int N) {

	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);

	int unsigned i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			lindblad_real[i + j * N] = 0.;
			lindblad_imag[i + j * N] = 0.;
		}
	}
	
	
	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, N);
}

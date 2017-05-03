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

// Stuff about MPI
#include <mpi.h>
#include <accel.h>  
#include <cudadevice.h>

//********************************************************************//
// defining global variables 

#define NSITES 6			 // number of excitonic sites to be modeled
#define dt 1.0				 // integration time step in fs
#define number_of_steps 10   // number of integration time steps
#define BLOCK_SIZE 16        // size for blocked matrix operations

#define TILE 16
#define BLOCK_SIZE 8

/*********************************************************************************/

extern void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda,
                                  double* w, double* work, int* lwork, int* info);

/*********************************************************************************/
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



// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
double _spectral_density(int j, double omega, double **params, int num_params) {
	double lambda, nu, omega_shift;
	double spec_dens; 
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	omega_shift = params[j][2];
	spec_dens  = nu * lambda * omega / (nu * nu + (omega - omega_shift) * (omega - omega_shift));
	spec_dens += nu * lambda * omega / (nu * nu + (omega + omega_shift) * (omega + omega_shift));
	return spec_dens;
}


double _spectral_density_derivative(int j, double omega, double **params, int num_params) {
// params ... [lambda_0, 1/nu_0, Omega_0, lambda_1, 1/nu_1, Omega_1, ...]
	double lambda, nu, omega_shift;
	lambda      = params[j][0];
	nu          = 1/params[j][1] * fs1_to_cm1;
	return 2 * lambda / nu;
}



double _phonon_statistics(double omega) {
	omega *= cm1_to_fs1;
	double result;
	result = 1. / (exp(HBAR * omega * BETA) - 1.);
	return result;
}



double _get_rate(int j, double omega, double **params, int num_params) {
	double result;
	if (omega < -1e-12) {
		result = 2 * PI * _spectral_density(j, -omega, params, num_params);
		result *= (_phonon_statistics(-omega) + 1) * cm1_to_fs1;
		return result;
	}
	if (omega > 1e-12) {
		result = 2 * PI * _spectral_density(j, omega, params, num_params);
		result *= _phonon_statistics(omega) * cm1_to_fs1;
		return result;
	}
//	printf("omega %.10f\n", omega);
	result = 2 * PI * KB * T * HBAR_INV;
	result *= _spectral_density_derivative(j, omega, params, num_params);
	return result;
}


void gen_test_spec_densities(double **params, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		params[i][0] = 35.;
		params[i][1] = 50.;
		params[i][2] = 0.;
	}
}



void get_rates(double ***gammas, double **params,  double *energies, int SIZE) {
	int unsigned j, M, N;
	double rate;
//	printf("getting rates\n");
	for (j = 0; j < SIZE; j++) {
		for (M = 0; M < SIZE; M++) {
			for (N = 0; N < SIZE; N++) {
//				printf("%d %d %d\n", j, M, N);
				if (j == 0) {
					gammas[M][N][j] = 0.;
//					printf("assigned zero\n");
				}
				if (j == SIZE - 1) {
//					printf("???\n");
					gammas[M][N][j] = 0.;
				}
				if ((j != 0) && (j != SIZE - 1)) {
//					printf("trying to get rate");
					rate = _get_rate(j - 1, energies[M] - energies[N], params, SIZE - 2);
//					printf("rate %d %d %d %.5f\n", j, M, N, rate);
					gammas[M][N][j] = rate;
				}
			}
		}
	}
}


void get_V_matrices(double ****V, double **eigVects, int N) {
	int unsigned i, j, k, l;
	double helper[N][N];
	// set everything to zero
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;
				V[k][j][i][1] = 0.;
				V[k][j][i][2] = 0.;
			}
		}
	}
	for (i = 1; i < N - 1; i++) {
		for (j = 0; j < N; j++) {
			for (k = 0; k < N; k++) {
				V[k][j][i][0] = 0.;
				V[k][j][i][1] = 0.;
				V[k][j][i][2] = 0.;
			}
		}
		// compute |m><m|
		V[i][i][i][1] = 1.;
		// rotate
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][1] * eigVects[l][k];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][1] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][1] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |0><m|
		V[0][i][i][0] = 1.;
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][0] * eigVects[l][k];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][0] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][0] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating

		// compute |RC><m|
		V[N - 1][i][i][2] = 1.;
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				helper[j][k] = 0.;
				for (l = 0; l < N; l++) {
					helper[j][k] += V[j][l][i][2] * eigVects[l][k];
				}
			}
		}
		for (k = 0; k < N; k++) {
			for (j = 0; j < N; j++) {
				V[j][k][i][2] = 0.;
				for (l = 0; l < N; l++) {
					V[j][k][i][2] += eigVects[l][j] * helper[l][k];
				}
			}
		} // done rotating
	}
}






void diagonalize(double **matrix, double *diagonal, int N) {
	int n = N, info, lwork;
	double wkopt;
	double* work;
	int LDA = N;

	int unsigned i, j;
	double *flat_matrix;
	flat_matrix = (double *) malloc(sizeof(double) * N*N);
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			flat_matrix[i + j * N] = matrix[i][j];
		}
	}

    // we need to get the optimal work load first
    lwork = -1;
    dsyev_("Vectors", "Upper", &n, flat_matrix, &LDA, diagonal, &wkopt, &lwork, &info);
    // then we perform the actual diagonalization
    lwork = (int)wkopt;
    work = (double*) malloc(lwork * sizeof(double));
    dsyev_("Vectors", "Upper", &n, flat_matrix, &LDA, diagonal, work, &lwork, &info);

    for (i = 0; i < N; i++) {
    	for (j = 0; j < N; j++) {
    		matrix[i][j] = flat_matrix[i + j * N];
    	}
    }
}


void gen_zero_matrix(double **A, int N) {
	int unsigned i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			A[i][j] = 0.;
		}
	}
}




// generate random hamiltonian
void gen_random_hamiltonian_real(double **H, int N) {
	int unsigned i, j;
	double number;
	gen_zero_matrix(H, N);
	for (i = 1; i < N - 1; i++) {
		for (j = 1; j < N - 1; j++) {
			if (i == j) 
				number = rand()%800 - 400;
			else
				number = rand()%200 - 100;
			H[i][j] = number;
			H[j][i] = number;
		}
	}
}



void gen_test_links(double *links_to_loss, double *links_to_target, int N) {
	int unsigned i;
	for (i = 0; i < N; i++) {
		links_to_loss[i] = 0.;
		links_to_target[i] = 0.;
	}
	for (i = 1; i < N - 1; i++) {
		links_to_loss[i] = 0.001;
		links_to_target[3] = 0.005;
	}
}



void transpose(double **matrix, int N) {
	int unsigned i, j;
	double element;
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i <= j) {
				element = matrix[i][j];
				matrix[i][j] = matrix[j][i];
				matrix[j][i] = element;
			}
		}
	}
}



void rotate(double **matrix, double **eigVects, int N) {
	int unsigned i, j, k;
	double helper[N][N];
//	printf("rotating\n");
	// A * eigVects
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			helper[i][j] = 0;
			for (k = 0; k < N; k++) {
				helper[i][j] += matrix[i][k] * eigVects[k][j];
			}
		}
	}
//	printf("second rotation\n");
	// eigVects^{-1} * A
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			matrix[i][j] = 0.;
			for (k = 0; k < N; k++) {
				matrix[i][j] += eigVects[k][i] * helper[k][j];
			}
		}
	}
//	printf("destroying helper\n");
//	free((void*) helper);
//	printf("ehlper destroyed\n");
}

// print matrix
void print_matrix_real(double **A, int N) {
	int unsigned i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; ++j) {
            printf("%.5f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double **alloc_2d_double(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

double ***alloc_3d_double(int rows, int cols, int dept) {
    double *data = (double *)malloc(rows*cols*dept*sizeof(double));
    double **array= (double **)malloc(rows*dept*sizeof(double*));
    double *** tensor = (double ***)malloc(dept*sizeof(double*));
    for (int i=0; i<dept; i++){
	tensor[i] = &(array[rows*i]);
    	for (int j=0; j<rows; j++){
        	tensor[i][j] = &(data[(rows*cols*i)+cols*j]);
	}
    }

    return tensor;
}

double ****alloc_4d_double(int rows, int cols, int dept,int side) {
    double *data = (double *)malloc(rows*cols*dept*side*sizeof(double));
    double **array= (double **)malloc(rows*dept*side*sizeof(double*));
    double *** tensor = (double ***)malloc(dept*side*sizeof(double*));
    double **** tensors = (double ***)malloc(side*sizeof(double*));
    for (int k=0; k<cols; k++){
	tensors[k] = &(tensor[dept*k]);
	for (int i=0; i<dept; i++){
		tensors[k][i] = &(array[(dept*cols*k)+cols*i]);
    		for (int j=0; j<side; j++){
        		tensors[k][i][j] = &(data[(rows*cols*side*k)+(rows*cols*i)+cols*j]);
		}
	}
    }

    return tensors;
}



int main(int argc, char *argv[]) {
	double tic, toc, total;
        // MPI stuff
        int myrank;
        int hostid;
	int ngpus,gpunum;
        
        int world_rank;
        int world_size;
        int number_amount;

        MPI_Status status;
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);


      	printf(" Initializing from node %d \n", world_rank);
        if (world_rank == 0) { // clock
    		/* Outputs */
   		FILE *tape;
    		tape = fopen("density.txt","a");
        
        }

	ngpus = acc_get_num_devices( acc_device_nvidia);
	if( ngpus ){
		gpunum = world_rank % ngpus;
		acc_set_device_num( gpunum, acc_device_nvidia );
	}
	else{
		 acc_set_device_type( acc_device_host );
	} 
	printf("\nNumgpus %d\n",ngpus);
	printf("\nNumgpus %d\n",acc_device_host);
	printf("\nNumgpus %d\n",gpunum);
	


        // Some variables for all nodes
	printf("# starting the propagation ...\n");

	// we start with initializing a number of helper integers and floats
	int unsigned i, j, k, l, SIZE;

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

	double **eigVects;
	double **rho_real;
	double **rho_imag;		// two dimensional arrays
	double **rho_real_host;
	double **rho_imag_host;		// two dimensional arrays
	
	eigVects = alloc_2d_double(SIZE,SIZE);
	rho_imag = alloc_2d_double(SIZE,SIZE);
	rho_real = alloc_2d_double(SIZE,SIZE);
	rho_imag_host = alloc_2d_double(SIZE,SIZE);
	rho_real_host = alloc_2d_double(SIZE,SIZE);

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... spectral density parameters (describe coupling to environment)
	// ... transition rates between excitonic states

	double **params;
	params = alloc_2d_double(SIZE-2,SIZE-2);

	double ***gammas;											// three dimensional array
	gammas = alloc_3d_double(SIZE,SIZE,SIZE);

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... storing transition matrices V
	// ... storing intermediate results

	double ****all_Vs;											// four dimensional array
	all_Vs = alloc_4d_double(SIZE,SIZE,SIZE,SIZE);

	double ***reduction_intermediates;							// three dimensional array
	reduction_intermediates = alloc_3d_double(SIZE,SIZE,SIZE);

	//------------------------------------------------------------------//

	// allocate space for ...
	// ... more intermediate results 
	// ... for computing the density matrix update
	// ... the density matrix update
	// ... the Runge Kutta steps

	// two dimensional arrays
	double **comm_real, **comm_imag, **lindblad_real, **lindblad_imag;
	comm_real = alloc_2d_double(SIZE,SIZE);
	comm_imag = alloc_2d_double(SIZE,SIZE);
	lindblad_real = alloc_2d_double(SIZE,SIZE);
	lindblad_imag = alloc_2d_double(SIZE,SIZE);

	// two dimensional arrays
	double **k1_real, **k2_real, **k3_real, **k4_real, **h1_real;
	k1_real = alloc_2d_double(SIZE,SIZE);
	k2_real = alloc_2d_double(SIZE,SIZE);
	k3_real = alloc_2d_double(SIZE,SIZE);
	k4_real = alloc_2d_double(SIZE,SIZE);
	h1_real = alloc_2d_double(SIZE,SIZE);

	double **k1_imag, **k2_imag, **k3_imag, **k4_imag, **h1_imag;
	k1_imag = alloc_2d_double(SIZE,SIZE);
	k2_imag = alloc_2d_double(SIZE,SIZE);
	k3_imag = alloc_2d_double(SIZE,SIZE);
	k4_imag = alloc_2d_double(SIZE,SIZE);
	h1_imag = alloc_2d_double(SIZE,SIZE);

	// three dimensional arrays
	double ***V, ***first, ***second, ***helper;
	V = alloc_3d_double(SIZE,SIZE,SIZE);
	first = alloc_3d_double(SIZE,SIZE,SIZE);
	second = alloc_3d_double(SIZE,SIZE,SIZE);
	helper = alloc_3d_double(SIZE,SIZE,SIZE);


	//------------------------------------------------------------------//

	// now we can set up our system

	// first, generate the hamiltonian
        if (world_rank == 0) { // clock
		gen_random_hamiltonian_real(eigVects, SIZE);
                for(i=1;i<world_size;i++){
       	        	MPI_Send(&(eigVects[0][0]),SIZE*SIZE, MPI_DOUBLE,i,10+i,  //Tag 0
               	      		MPI_COMM_WORLD); 
        	}
	}
	else{
      		printf(" I am copying %d \n", world_rank);
		MPI_Recv(&(eigVects[0][0]),SIZE*SIZE, MPI_DOUBLE, 0, 10+world_rank, MPI_COMM_WORLD,
 			&status);

    		MPI_Get_count(&status, MPI_DOUBLE, &number_amount);
                printf("1 received %d numbers from 0. Message source = %d, "
                       "tag = %d /n",
                       number_amount, status.MPI_SOURCE, status.MPI_TAG);
        }
        print_matrix_real(eigVects,SIZE);

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
	for (i = 0; i < SIZE; i++) {
		for (j = 0; j < SIZE; j++) {
			lindblad_real[i][j] = 0.;
			lindblad_imag[i][j] = 0.;
			for (k = 0; k < SIZE; k++) {
				reduction_intermediates[i][j][k] = 0.0;

			}
		}
	}


	// start propagation //
	int unsigned ii, jj;
	int unsigned index, jndex, kndex;
	int unsigned step;
	tic = clock();
	int var, value;
	value = 29;
	int unsigned m, M, N;
	int unsigned kk;
	double rate, sum;	
	double iupper, jupper;


	// copy all data to the GPU
	// ... note: once all the data is on the GPU, only the density matrix needs to be communicated between GPU and CPU
	#pragma acc enter data copyin(hamiltonian[0:SIZE], gammas[0:SIZE][0:SIZE][0:SIZE], all_Vs[0:SIZE][0:SIZE][0:SIZE][0:3], links_to_target[0:SIZE], links_to_loss[0:SIZE])  
	#pragma acc enter data copyin(V[0:SIZE][0:SIZE][0:SIZE], first[0:SIZE][0:SIZE][0:SIZE], second[0:SIZE][0:SIZE][0:SIZE], helper[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc enter data copyin(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], k2_real[0:SIZE][0:SIZE], k2_imag[0:SIZE][0:SIZE], k3_real[0:SIZE][0:SIZE], k3_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
	#pragma acc enter data copyin(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE], lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE])
	#pragma acc enter data copyin(reduction_intermediates[0:SIZE][0:SIZE][0:SIZE])
	#pragma acc enter data copyin(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])
	for (step = 0; step < number_of_steps; step++) {


		get_density_update(rho_real, rho_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, &lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);
		#pragma acc update host(lindblad_real[0:SIZE][0:SIZE])
		MPI_Allreduce(MPI_IN_PLACE, &(lindblad_real[0][0]),SIZE*SIZE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		#pragma acc update device(lindblad_real[0:SIZE][0:SIZE])
        	
		// ... compute k1 in blocked matrix operations
		#pragma acc kernels present(comm_real[0:SIZE][0:SIZE], comm_imag[0:SIZE][0:SIZE]) copyin(lindblad_real[0:SIZE][0:SIZE], lindblad_imag[0:SIZE][0:SIZE]) present(rho_real[0:SIZE][0:SIZE], rho_imag[0:SIZE][0:SIZE])     present(k1_real[0:SIZE][0:SIZE], k1_imag[0:SIZE][0:SIZE], h1_real[0:SIZE][0:SIZE], h1_imag[0:SIZE][0:SIZE])
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

		get_density_update(h1_real, h1_imag, hamiltonian, comm_real, comm_imag, gammas, eigVects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, V, first, second, helper, reduction_intermediates, SIZE);
		#pragma acc update host(lindblad_real[0:SIZE][0:SIZE])
		MPI_Allreduce(MPI_IN_PLACE, &(lindblad_real[0][0]),SIZE*SIZE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		#pragma acc update device(lindblad_real[0:SIZE][0:SIZE])
                
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

		#pragma acc update host(rho_real[0:SIZE][0:SIZE])


		transpose(eigVects, SIZE);
		rotate(rho_real, eigVects, SIZE);
		rotate(rho_imag, eigVects, SIZE);
		transpose(eigVects, SIZE);


	}
	total = 0;
	printf("%d ", step);
        for (i = 0; i < SIZE; i++) {
            printf("%.10f ", rho_real[i][i]);
		total += rho_real[i][i];
        }


    // Analyze time elapsed
    //double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    //printf("# Time Elapsed: %f seconds\n\n", time_spent);
}


#include <stdio.h>
#include <stdlib.h>

#include "headers.h"

/* methods for computing the Liouville operator on the density matrix
 * ... all methods assume operators to be in the exciton basis, unless otherwise specified
 */

#define cm1_to_fs1 1. / 33356.40952
#define fs1_to_cm1 33356.40952
#define eV_to_cm1  8065.54429
#define s_to_fs    1.e15

#define HBAR 6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV 1. / (6.582119514e-16*eV_to_cm1*s_to_fs)	// FIXME: Units might be off!

/********************************************************************/

 void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N) {
 	// FIXME
 	// the hamiltonian is a real valued diagonal matrix
 	// ... but we don't care for now and use simple matrix multiplication

 	double *h_real, *h_imag; 
 	h_real = (double *) malloc(sizeof(double) * (N * N));
 	h_imag = (double *) malloc(sizeof(double) * (N * N));
 	gen_zero_matrix_real(h_real, N);
 	gen_zero_matrix_real(h_imag, N);

 	double *help1_real, *help1_imag, *help2_real, *help2_imag;
 	help1_real = (double *) malloc(sizeof(double) * (N * N));
	help2_real = (double *) malloc(sizeof(double) * (N * N));
 	help1_imag = (double *) malloc(sizeof(double) * (N * N));
	help2_imag = (double *) malloc(sizeof(double) * (N * N));

	// get the first part of the commutator
	int unsigned i;
	for (i = 0; i < N; i++) {
		h_real[i + i * N] = hamiltonian[i];
	}

	matrix_mul_complexified(h_real, h_imag, rho_real, rho_imag, help1_real, help1_imag, N);
	// get the second part
	matrix_mul_complexified(rho_real, rho_imag, h_real, h_imag, help2_real, help2_imag, N);

	// combine the two parts
	matrix_sub_real(help1_imag, help2_imag, comm_real, N);
	matrix_sub_real(help2_real, help1_real, comm_imag, N);


//	matrix_mul_complexified(hamiltonian, h_imag, rho_real, rho_imag, help1_real, help1_imag, N);
	// get the second part of the commutator
//	matrix_mul_complexified(rho_real, rho_imag, hamiltonian, h_imag, help2_real, help2_imag, N);
	// combine the two parts
//	matrix_add_real(help1_imag, help2_imag, comm_real, N);
//	matrix_add_real(help2_real, help1_real, comm_imag, N);
	// don't forget about hbar

	matrix_mul_scalar(comm_real, HBAR_INV, N);
	matrix_mul_scalar(comm_imag, HBAR_INV, N);

	free((void*) h_real);
	free((void*) h_imag);
	free((void*) help1_real);
	free((void*) help1_imag);
	free((void*) help2_real);
	free((void*) help2_imag);

 }



void lindblad_operator(double *rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, int SIZE) {
	int unsigned m, M, N;
	double rate;

	double *V;
	V = (double *) malloc(sizeof(double) * SIZE * SIZE);
	double *V_dagg;
	V_dagg = (double *) malloc(sizeof(double) * SIZE * SIZE);

	double *first_real, *second_real, *third_real, *helper;//, *first_imag, *second_imag, *third_imag;
	first_real  = (double *) malloc(sizeof(double) * SIZE * SIZE);
	second_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	third_real  = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper      = (double *) malloc(sizeof(double) * SIZE * SIZE);


	// decay of excitons into other excitonic states
	for (m = 1; m < SIZE - 1; m++) {
		for (M = 0; M < SIZE; M++) {
			for (N = 0; N < SIZE; N++) {
				rate = gammas[M + N * SIZE + m * SIZE * SIZE];
//				rate = (double)(rate);

				get_V(V, eigVects, m, m, SIZE);
				get_V(V_dagg, eigVects, m, m, SIZE);
				transpose(V_dagg, SIZE);

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

				matrix_add_real(lindblad_real, first_real, lindblad_real, SIZE);
		
			}
		}
	}


	// decay of excitons into the ground (loss) state
	for (m = 0; m < SIZE; m++) {
		rate = links_to_loss[m];

		get_V(V, eigVects, 0, m, SIZE);
		get_V(V_dagg, eigVects, 0, m, SIZE);
		transpose(V_dagg, SIZE);

		matrix_mul_real(V, rho_real, helper, SIZE);
		matrix_mul_real(helper, V_dagg, first_real, SIZE);

		matrix_mul_real(rho_real, V_dagg, helper, SIZE);
		matrix_mul_real(helper, V, second_real, SIZE);
		matrix_mul_scalar(second_real, 0.5, SIZE);

		matrix_mul_real(V_dagg, V, helper, SIZE);
		matrix_mul_real(helper, rho_real, third_real, SIZE);
		matrix_mul_scalar(third_real, 0.5, SIZE);	

		matrix_sub_real(first_real, second_real, first_real, SIZE);
		matrix_sub_real(first_real, third_real, first_real, SIZE);
		matrix_mul_scalar(first_real, rate, SIZE);

		matrix_add_real(lindblad_real, first_real, lindblad_real, SIZE);		
	}

//	print_matrix_real(rho_real, SIZE);

	// decay of excitons into the target (reaction center) state
	for (m = 0; m < SIZE; m++) {
		rate = links_to_target[m];

		get_V(V, eigVects, SIZE - 1, m, SIZE);
		get_V(V_dagg, eigVects, SIZE - 1, m, SIZE);
		transpose(V_dagg, SIZE);

		matrix_mul_real(V, rho_real, helper, SIZE);
		matrix_mul_real(helper, V_dagg, first_real, SIZE);

		matrix_mul_real(V_dagg, V, helper, SIZE);
		matrix_mul_real(helper, rho_real, second_real, SIZE);
		matrix_mul_scalar(second_real, 0.5, SIZE);

		matrix_mul_real(rho_real, V_dagg, helper, SIZE);
		matrix_mul_real(helper, V, third_real, SIZE);
		matrix_mul_scalar(third_real, 0.5, SIZE);	

		matrix_sub_real(first_real, second_real, first_real, SIZE);
		matrix_sub_real(first_real, third_real, first_real, SIZE);
		matrix_mul_scalar(first_real, rate, SIZE);

		matrix_add_real(lindblad_real, first_real, lindblad_real, SIZE);		
	}


	free((void*) V);
	free((void*) V_dagg);
	free((void*) first_real);
	free((void*) second_real);
	free((void*) third_real);
	free((void*) helper);

}




void get_density_update(double *rho_real, double *rho_imag, double *energies, double *comm_real, double *comm_imag, 
	                    double *gammas, double *eigvects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, int N) {

	hamiltonian_commutator(rho_real, rho_imag, energies, comm_real, comm_imag, N);
	gen_zero_matrix_complex(lindblad_real, lindblad_imag, N);
	lindblad_operator(rho_real, rho_imag, gammas, eigvects, lindblad_real, lindblad_imag, links_to_loss, links_to_target, N);

}
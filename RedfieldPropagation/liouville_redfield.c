
#include <stdio.h>
#include <stdlib.h>

#include "headers.h"

/* methods for computing the Liouville operator on the density matrix
 * ... all methods assume operators to be in the exciton basis, unless otherwise specified
 */

#define cm1_to_fs1 1. / 33356.40952
#define fs1_to_cm1 33356.40952
#define eV_to_cm1  8065.54429
#define s_to_fs    10.e15

#define HBAR 6.582119514e-16*eV_to_cm1*s_to_fs
#define HBAR_INV 10. / (6.582119514e-16*eV_to_cm1*s_to_fs)	// FIXME: Units might be off!

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

 }



void lindblad_operator(double *rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, int SIZE) {
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
//	first_imag  = (double *) malloc(sizeof(double) * SIZE * SIZE);
//	second_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
//	third_imag  = (double *) malloc(sizeof(double) * SIZE * SIZE);

	for (m = 1; m < SIZE - 1; m++) {
		for (M = 1; M < SIZE - 1; M++) {
			for (N = 1; N < SIZE - 1; N++) {
				rate = gammas[M + N * SIZE + m * SIZE * SIZE];
				rate = (double)(rate);
//				printf("%.10f\n", rate);
//				exit(1);

				get_V(V, eigVects, m, m, SIZE);
				get_V(V_dagg, eigVects, m, m, SIZE);
				transpose(V_dagg, SIZE);

//				print_matrix_real(V, SIZE);
//				print_matrix_real(V_dagg, SIZE);

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
//				print_matrix_real(first_real, SIZE);
//				printf("LINDBLAD\n");
//				print_matrix_real(lindblad_real, SIZE);
			}
		}
	}
}
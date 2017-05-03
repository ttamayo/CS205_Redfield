

// external modules
#include <stdio.h>
#include <stdlib.h>

// internal modules
#include "headers.h"

/**********************************************************************/

#define NSITES 14
#define dt 1.0
#define number_of_steps 1

/**********************************************************************/


int main(void) {
	printf("# starting\n");
	int unsigned i, j, k;
	double tic, toc;
	int SIZE;
	SIZE = NSITES + 2;
	printf("# SIZE %d\n", SIZE);

	double *A, *D, *gammas, *params, *links_to_loss, *links_to_target;
	A               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	D               = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gammas          = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE);
	params          = (double *) malloc(sizeof(double) * 3 * (SIZE - 2));
	links_to_loss   = (double *) malloc(sizeof(double) * SIZE);
	links_to_target = (double *) malloc(sizeof(double) * SIZE);


	gen_random_hamiltonian_real(A, SIZE);
	gen_test_spec_densities(params, NSITES);
	gen_test_links(links_to_loss, links_to_target, SIZE);

	double *rho_real, *rho_imag, *helper_matrix;
	rho_real      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	rho_imag      = (double *) malloc(sizeof(double) * SIZE * SIZE);
	helper_matrix = (double *) malloc(sizeof(double) * SIZE * SIZE);
	gen_zero_matrix_complex(rho_real, rho_imag, SIZE);
	rho_real[1 + 1 * SIZE] = 1.;

	double *all_Vs;
	all_Vs = (double *) malloc(sizeof(double) * SIZE*SIZE*SIZE*3);
	// diagonalize hamiltonian
	diagonalize(A, D, SIZE);
	

	// get rates in site basis
	get_V_matrices(all_Vs, A, SIZE);	
	get_rates(gammas, params, D, NSITES, SIZE);


	// rotate rho into exciton basis
	rotate(rho_real, A, SIZE);
	rotate(rho_imag, A, SIZE);

	double *comm_real, *comm_imag;
	comm_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	comm_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	double *lindblad_real, *lindblad_imag;
	lindblad_real = (double *) malloc(sizeof(double) * SIZE * SIZE);
	lindblad_imag = (double *) malloc(sizeof(double) * SIZE * SIZE);
	

 
	get_density_update(rho_real, rho_imag, D, comm_real, comm_imag, gammas, A, lindblad_real, lindblad_imag, links_to_loss, links_to_target, all_Vs, SIZE);


        for (i = 0; i < SIZE; i++) {
       		printf("%.10f ", rho_real[i + i * SIZE]);
        }
        printf("\n");


	free((void*) comm_real);
	free((void*) comm_imag);
	free((void*) lindblad_real);
	free((void*) lindblad_imag);



	return 0;
}

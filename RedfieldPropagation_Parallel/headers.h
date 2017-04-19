
<<<<<<< HEAD
#ifndef FUNC_LIOUVILLE_INCLUDED
	#define FUNC_LIOUVILLE_INCLUDED

	void get_density_update(double *rho_real, double *rho_imag, double *energies, double *comm_real, double *comm_imag, 
	                    double *gammas, double *eigvects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int N);

	#pragma acc routine gang
	void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N);
	#pragma acc routine gang
	void lindblad_operator(double * restrict rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, double *all_Vs, int SIZE);
=======
<<<<<<< HEAD
#ifndef FUNC_LIOUVILLE_INCLUDED
	#define FUNC_LIOUVILLE_INCLUDED

	#pragma acc routine gang
	void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N);
	#pragma acc routine gang
	void lindblad_operator(double *rho_real, double *rho_imag, double *gammas, double *eigVects, double *lindblad_real, double *lindblad_imag, double *links_to_loss, double *links_to_target, int SIZE);
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
#endif

/******************************************************************************************/

#ifndef FUNC_RATES_INCLUDED
	#define FUNC_RATES_INCLUDED
	
	void get_rates(double *gammas, double *params, double *energies, int num_params, int Nsites2);
<<<<<<< HEAD
//	#pragma acc routine worker
	void get_V(double *V, double *eigvects, int i, int k, int N);
	#pragma acc routine gang
	void get_V_matrices(double *V, double *eigvects, int N);
=======
	#pragma acc routine worker
	void get_V(double *V, double *eigvects, int i, int k, int N);
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
#endif

/******************************************************************************************/
#ifndef FUNC_MMO_INCLUDED
	#define FUNC_MMO_INCLUDED

<<<<<<< HEAD
//	#pragma acc routine worker
	void matrix_add_real(double *A, double *B, double *C, int N);

//	#pragma acc routine worker
	void matrix_sub_real(double *A, double *B, double *C, int N);

//	#pragma acc routine worker
	void matrix_mul_real(double *A, double *B, double *C, int N);

//	#pragma acc routine worker
	void matrix_mul_scalar(double *A, double scalar, int N);

	void matrix_add_complex(double *A_real, double *A_imag, double *B_real, double *B_imag, double *C_real, double *C_imag, int N);
=======
	#pragma acc routine worker
	void matrix_add_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_sub_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_mul_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_mul_scalar(double *A, double scalar, int N);
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
#endif

/******************************************************************************************/


#ifndef FUNC_MAT_GEN_INCLUDED
	#define FUNC_MAT_GEN_INCLUDED

	void gen_identity_matrix(double *A, int N);
<<<<<<< HEAD
//	#pragma acc routine worker
=======
	#pragma acc routine worker
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	void gen_zero_matrix_real(double *A, int N);
	void gen_one_matrix_real(double *A, int N);
	void gen_random_matrix_real(double *mat, int N);
	void gen_random_hamiltonian_real(double *H, int N);
	void gen_test_hamiltonian(double *A);
	void gen_test_links(double *links_to_loss, double *links_to_target, int N);
	void gen_test_spec_densities(double *params, int N);
	void gen_identity_complex(double *A_real, double *A_imag, int N);
<<<<<<< HEAD
//	#pragma acc routine worker
=======
	#pragma acc routine worker
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa
	void gen_zero_matrix_complex(double *A_real, double *A_imag, int N);


#endif


#ifndef FUNC_SMO_INCLUDED
	#define FUNC_SMO_INCLUDED
	// from single_matrix_operations.c
	void diagonalize(double *A, double *D, int N);
<<<<<<< HEAD
//	#pragma acc routine worker
	void rotate(double *A, double *eigvect, int N);
//	#pragma acc routine worker
	void transpose(double *A, int N);
#endif



#ifndef FUNC_UTIL_INCLUDED
	#define FUNC_UTIL_INCLUDED
//	#pragma acc routine
	void print_matrix_real(double *A, int N);
#endif


=======
	#pragma acc routine worker
	void rotate(double *A, double *eigvect, int N);
	#pragma acc routine worker
	void transpose(double *A, int N);
#endif
=======
// from liouville_redfield.c

void hamiltonian_commutator(double *rho_real, double *rho_imag, double *hamiltonian, double *comm_real, double *comm_imag, int N);

// from matrix_generators.c
void gen_identity_matrix(double *A, int N);
#pragma acc routine
void gen_zero_matrix_real(double *A, int N);
void gen_one_matrix_real(double *A, int N);
void gen_random_matrix_real(double *mat, int N);
void gen_random_hamiltonian_real(double *H, int N);
void gen_test_hamiltonian(double *A);
void gen_test_links(double *links_to_loss, double *links_to_target, int N);
void gen_test_spec_densities(double *params, int N);
void gen_identity_complex(double *A_real, double *A_imag, int N);
void gen_zero_matrix_complex(double *A_real, double *A_imag, int N);


// from multiple_matrix_operations.c
#pragma acc routine
void matrix_add_real(double *A, double *B, double *C, int N);
#pragma acc routine 
void matrix_sub_real(double *A, double *B, double *C, int N);
#pragma acc routine
void matrix_mul_real(double *A, double *B, double *C, int N);
#pragma acc routine
void matrix_mul_scalar(double *A, double scalar, int N);

// from rates.c
void get_rates(double *gammas, double *params, double *energies, int num_params, int Nsites2);
#pragma acc routine
void get_V(double *V, double *eigvects, int i, int k, int N);

// from single_matrix_operations.c
void diagonalize(double *A, double *D, int N);
#pragma acc routine
void rotate(double *A, double *eigvect, int N);
#pragma acc routine
void transpose(double *A, int N);


// from utilities.c
void print_matrix_real(double *A, int N);
>>>>>>> 8a8db8f89679f383fb17ba16f8a581ae653ec48e
>>>>>>> 4dd0983a47a7cf0708d39e1c642cc42d3cd430fa

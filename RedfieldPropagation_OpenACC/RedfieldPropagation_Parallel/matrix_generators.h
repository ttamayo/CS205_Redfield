
#ifdef FUNCTIONS_H_INCLUDED 
	#define FUNCTIONS_H_INCLUDED
	
	void gen_identity_matrix(double *A, int N);
	#pragma acc routine vector
	void gen_zero_matrix_real(double *A, int N);
	void gen_one_matrix_real(double *A, int N);
	void gen_random_matrix_real(double *mat, int N);
	void gen_random_hamiltonian_real(double *H, int N);
	void gen_test_hamiltonian(double *A);
	void gen_test_links(double *links_to_loss, double *links_to_target, int N);
	void gen_test_spec_densities(double *params, int N);
	void gen_identity_complex(double *A_real, double *A_imag, int N);
	#pragma acc routine vector
	void gen_zero_matrix_complex(double *A_real, double *A_imag, int N);

#endif

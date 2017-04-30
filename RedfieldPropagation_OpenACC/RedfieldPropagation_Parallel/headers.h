


#ifndef FUNC_LIOUVILLE_INCLUDED
	#define FUNC_LIOUVILLE_INCLUDED
	void hamiltonian_commutator(double **rho_real, double **rho_imag, double *hamiltonian, double **comm_real, double **comm_imag, int N);
	void lindblad_operator(double **rho_real, double **rho_imag, double ***gammas, double **eigVects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ** V, double **first_real, double **second_real, double **helper, double ***reduction_intermediates,  int SIZE);
	void get_density_update(double **rho_real, double **rho_imag, double *energies, double **comm_real, double **comm_imag, double ***gammas, double **eigvects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double **V, double **first_real, double **second_real, double **helper, double ***reduction_intermediates, int N); 
#endif



#ifndef FUNC_UTIL_INCLUDED
	#define FUNC_UTIL_INCLUDED
	void gen_zero_matrix(double **A, int N); 
	void gen_random_hamiltonian_real(double **H, int N); 
	void gen_test_links(double *links_to_loss, double *links_to_target, int N);
	void gen_test_spec_densities(double **params, int N);
	void print_matrix_real(double **A, int N);
	void rotate(double **matrix, double **eigVects, int N);
	void transpose(double **matrix, int N);
#endif



#ifndef FUNC_RATES_INCLUDED
	#define FUNC_RATES_INCLUDED
	void get_rates(double ***gammas, double **params,  double *energies, int SIZE);
	void get_V_matrices(double ****V, double **eigVects, int N);
#endif

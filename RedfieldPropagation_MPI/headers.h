


#ifndef FUNC_LIOUVILLE_INCLUDED
	#define FUNC_LIOUVILLE_INCLUDED
	void hamiltonian_commutator(double **rho_real, double **rho_imag, double *hamiltonian, double **comm_real, double **comm_imag, int N);
	void lindblad_operator(double **rho_real, double **rho_imag, double ***gammas, double **eigVects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates,  int SIZE);
	void get_density_update(double **rho_real, double **rho_imag, double *energies, double **comm_real, double **comm_imag, double ***gammas, double **eigvects, double **lindblad_real, double **lindblad_imag, double *links_to_loss, double *links_to_target, double ****all_Vs, double ***V, double ***first, double ***second, double ***helper, double ***reduction_intermediates, int N); 
#endif





#ifndef FUNC_RATES_INCLUDED
	#define FUNC_RATES_INCLUDED
	void get_rates(double ***gammas, double **params,  double *energies, int SIZE);
	void get_V_matrices(double ****V, double **eigVects, int N);
#endif

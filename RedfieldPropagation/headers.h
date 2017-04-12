

// from matrix_generators.c
void gen_identity_matrix(double *A, int N);
void gen_zero_matrix_real(double *A, int N);
void gen_one_matrix_real(double *A, int N);
void gen_test_hamiltonian(double *A);

void gen_identity_complex(double *A_real, double *A_imag, int N);


// from multiple_matrix_operations.c
void matrix_mul_real(double *A, double *B, double *C, int N);
void matrix_mul_scalar(double *A, double scalar, int N);

// from rates.c
void get_rates(double *gammas, double *params, double *energies, int num_params, int Nsites2);
void get_V(double *V, double *eigvects, int i, int k, int N);

// from single_matrix_operations.c
void diagonalize(double *A, double *D, int N);
void rotate(double *A, double *eigvect, int N);
void transpose(double *A, int N);


// from utilities.c
void print_matrix_real(double *A, int N);
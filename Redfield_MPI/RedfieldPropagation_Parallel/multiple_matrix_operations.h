
#ifdef FUNC_MMO_INCLUDED
	#define FUNC_MMO_INCLUDED

	#pragma acc routine worker
	void matrix_add_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_sub_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_mul_real(double *A, double *B, double *C, int N);

	#pragma acc routine worker
	void matrix_mul_scalar(double *A, double scalar, int N);

#endif

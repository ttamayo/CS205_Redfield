
#ifdef FUNC_MMO_INCLUDED
	#define FUNC_MMO_INCLUDED
	void matrix_add_real(double *A, double *B, double *C, int N);
	void matrix_sub_real(double *A, double *B, double *C, int N);
	void matrix_mul_real(double *A, double *B, double *C, int N);
	void matrix_mul_scalar(double *A, double scalar, int N);
#endif

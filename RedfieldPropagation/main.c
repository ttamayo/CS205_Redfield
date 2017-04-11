

// external modules
#include <stdio.h>
#include <stdlib.h>

// internal modules
#include "headers.h"

/**********************************************************************/


/**********************************************************************/


int main(void) {

	double *A;
	A = (double *) malloc(sizeof(double) * 4);

	gen_identity_real(A, 4);
	print_matrix_real(A, 4);

	return 0;

	// FIXME
	// not sure why this should be after the return command
	// putting this line before the return raises an error 
	free((void*) A);

}
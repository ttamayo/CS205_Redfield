

// external modules

#include <stdio.h>
#include <stdlib.h>

// internal modules

#include "matrix_generators.h"
#include "utilities.h"

/**********************************************************************/


/**********************************************************************/


int main() {

	double A[4];
	gen_identity_real(A, 4);
	print_matrix_real(A, 4);
	return 0;
}
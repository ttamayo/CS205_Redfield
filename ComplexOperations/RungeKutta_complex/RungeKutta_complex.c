/*
 * C program for implementation of Runge Kutta 4th order integration
 * Complex implementation using standard library.
 *
 * For our problem, \rho = y
 *
 * author: hannahsim
 *
 */ 

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// matrix size: SIZE * SIZE
#define SIZE 64

// timestep (in fs)
#define TIMESTEP 0.2


// Rate function, dy/dx (testing for now...)
double complex rate(double complex x, double complex y) {
    return -2.*y+x+4.;
}


// For testing purposes (initializing matrices)
void get_random_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            // each component between 0 and 9
            mat[i + j * N] = rand()%10 + rand()%10*I;
}

void get_zero_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 0. + 0.*I;
}

void get_one_matrix(double complex *mat, int N) {
    int unsigned i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            mat[i + j * N] = 1. + 0.*I;
}




// Helper function: individual element Runge-Kutta 4th order
// returns y(x+h) i.e. y after taking a single time step h
double complex RK4_element(double complex x, double complex y, double dx) {     
    double complex k1 = dx * rate(x,y),    
                   k2 = dx * rate(x + dx / 2. , y + k1 / 2.),
                   k3 = dx * rate(x + dx / 2. , y + k2 / 2.),
                   k4 = dx * rate(x + dx, y + k3); 
    
    return y + (k1 + 2.*k2 + 2.*k3 + k4) / 6.;
} 




// 4th order Runge Kutta for matrix over (n_end-n_init) number of steps with each stepsize = dx.
void RK4_matrix(double complex *xmat, double complex *ymat, int n_init, int n_end, double dx, int N) {
    
    int unsigned i,j,t;
    
    for (t=n_init; t<n_end; ++t) { 
        for (i=0;i<N;++i) {
            for (j=0;j<N;++j) {
                ymat[i+j*N] = RK4_element(xmat[i+j*N],ymat[i+j*N],dx);
            }
        }

/*        
        // show matrix after every 20 iterations
        if (t%20 == 0) {
            printf("Iteration %i\n",t);
            for(int k = 0; k < SIZE*SIZE; k++) {
                printf("%f + %fi ", creal(ymat[k]),cimag(ymat[k]));
                if(((k + 1) % SIZE) == 0)
                printf("\n");
            }
        }
*/

    }
}



// PROGRAM MAIN
int main() {

    // Allocate matrices
    double complex *ymat;
    double complex *xmat;
    ymat = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));
    xmat = (double complex *) malloc(sizeof(double complex) * (SIZE * SIZE));

    // Initialize matrices
    get_random_matrix(ymat, SIZE);  
    get_zero_matrix(xmat, SIZE);

    // Setup for timing
    clock_t tic, toc;
    tic = clock();

    int n_0 = 0;
    int n_f = 10; // number of steps (stepsize = TIMESTEP or "h")
 
    printf("\n\n");
    printf("****************************\n");
    printf("RUNGE-KUTTA4 IMPLEMENTATION:\n");
    printf("****************************\n");
 
    RK4_matrix(xmat,ymat,n_0,n_f,TIMESTEP,SIZE);
    
    toc = clock();

    // Analyze time elapsed
    double time_spent = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\nRESULTS:\n");
    printf("--------\n");
    printf("Time Elapsed: %f seconds\n\n", time_spent);


/*
    // Show final matrix (only useful for testing small enough matrices)
    printf("\nMatrix y_final\n");
    for(int i = 0; i < SIZE*SIZE; i++) {
        printf("%f + %fi ", creal(ymat[i]),cimag(ymat[i]));
        if(((i + 1) % SIZE) == 0)
        printf("\n");
    }
    printf("\n\n");    

*/

    return 0;
}

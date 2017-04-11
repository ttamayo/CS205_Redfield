/*!
 * This is a 
 * 
 *
 *
 * Versions:
 *    Teresa Tamayo - 03/30/2107
 */
#include <stdio.h>
#include <math.h>       
#include <time.h>
#include <string.h>

/* Size of the hamiltonian */
#define N 6
#define LDA N
#define LDVL N
#define LDVR N


/* Propagations */
#define DT 0.25
#define STEPS 200

/* Spectral density */
#define NUMBER 1
#define lambda_k[NUMBER] = {35.5f};
#define nu_k[1] = {35.5f};

/*
 * Physical conversions
 */
#define S_TO_FS 10.0e15  
#define EV_TO_CM1 8065.54429
#define FS1_TO_CM1 33356.40952
#define CM1_TO_FS1  1.0/33356.40952

/*
 * Physical constants
 */
#define HBAR  6.582119514e-16*EV_TO_CM1*S_TO_FS    //                # cm-1 * fs
#define T 300.                                    //                        # K
#define KB  8.6173303e-5*EV_TO_CM1                //                 # cm1 / K
#define BETA 1.0/(KB * T)                        //           # 1 / cm-1

/* Complex datatype */
struct _fcomplex {float re, im; };
typedef struct _fcomplex fcomplex;

/* CGEEV prototype */
extern void cgeev( char* jobvl, char* jobvr, int* n, fcomplex* a,
                int* lda, fcomplex* w, fcomplex* vl, int* ldvl, fcomplex* vr, int* ldvr,
                fcomplex* work, int* lwork, float* rwork, int* info );

void *memcpy(void *dest, const void *src, size_t n);
/* */
#define IDX(i,j,num)(i+j*num)

//#include <random>
#ifndef _LIMITS_H___
#define _LIMITS_H___

/* Number of bits in a `char'.  */
#undef CHAR_BIT
#define CHAR_BIT __CHAR_BIT__

/* Maximum length of a multibyte character.  */
#ifndef MB_LEN_MAX
#define MB_LEN_MAX 1
#endif

/* Minimum and maximum values a `signed char' can hold.  */
#undef SCHAR_MIN
#define SCHAR_MIN (-SCHAR_MAX - 1)
#undef SCHAR_MAX
#define SCHAR_MAX __SCHAR_MAX__

/* Maximum value an `unsigned char' can hold.  (Minimum is 0).  */
#undef UCHAR_MAX
#if __SCHAR_MAX__ == __INT_MAX__
# define UCHAR_MAX (SCHAR_MAX * 2U + 1U)
#else
# define UCHAR_MAX (SCHAR_MAX * 2 + 1)
#endif
#endif /* _LIMITS_H___ */


/*
 * Priting matrix
 */
void print_matrix( char* desc, int m, int n, fcomplex* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
                printf( "\n" );
        }
}

/* 
 * Function to get random seeds
 */

unsigned time_seed()
{
    time_t now = time(0);
    unsigned char *p = (unsigned char *)&now;
    unsigned seed = 0;
    size_t i;

    for (i = 0; i < sizeof now; i++)
    {
        seed = seed * (UCHAR_MAX) + p[i];
    }

    return seed;
}

/*
 * This function builds H
 */
void get_initialrho(fcomplex *p,int m)
{    
  int i,j;
  srand(time_seed());
  for(i=0; i<m; i++){
       for(j=0; j<m; j++){
           p[IDX(i,j,m)].re = 0.0;
           p[IDX(i,j,m)].im = 0.0;
       }
  } 
  p[1,1].re = 1.0;
}

/* 
 * Diagonal H
 */
void get_diaghamiltonian(fcomplex *h, fcomplex *eigen, int m)
{
  int i,j;
  for(i=0; i<m; i++){
       for(j=0; j<m; j++){
           h[IDX(i,j,m)].re = 0.0f;
           h[IDX(i,j,m)].im = 0.0f;
       }
       h[IDX(i,j,m)].re = eigen[0].re;
       h[IDX(i,j,m)].im = eigen[1].im;
  } 
}


/*
 * This function builds H
 */
void get_hamiltonian(fcomplex *h,int m)
{    
  int i,j;
  srand(time_seed());
  for(i=0; i<m; i++){
       for(j=0; j<m; j++){
           h[IDX(i,j,m)].re = 0.0f;
           h[IDX(i,j,m)].im = 0.0f;
       }
  } 
  // Example hamiltonian
  h[IDX(1,4,m)].re = 416.40098494f;
  h[IDX(2,4,m)].re = 193.24458563f;
  h[IDX(3,4,m)].re = 1228.67310624f;
  h[IDX(4,4,m)].re = -386.50195718f;

  h[IDX(1,1,m)].re = -29.79868321f;
  h[IDX(1,2,m)].re = 554.42775405f;
  h[IDX(1,3,m)].re = 584.77457258f;
  h[IDX(2,2,m)].re = -106.0223773f;
  h[IDX(2,3,m)].re = 325.25929639f;
  h[IDX(3,3,m)].re = -582.96766952f;
  h[IDX(4,3,m)].re = 1228.67310624f;
  h[IDX(2,1,m)].re = 554.42775405f;
  h[IDX(3,1,m)].re = 584.77457258f;
  h[IDX(4,1,m)].re = 416.40098494f;
  h[IDX(3,2,m)].re = 325.25929639f;
  h[IDX(4,2,m)].re = 193.24458563f;

}

/* 
 * Matrix multiplication
 */

void blas_matrix_multiplication(fcomplex *a, fcomplex *b, fcomplex *c, fcomplex alpha, fcomplex beta, int n)
{
     int k=N;
     int m=N;
     int lda = N;
     int ldb = N;
     int ldc = N;
     cgemm("N","N",&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc);
     return;
}


/*
 * Inverse of a matrix
 */

void inversematrix(fcomplex *a, int n)
{   
    int lda = N;
    int ipiv[N];
    int info;
    int lwork = N*N;
    fcomplex works[N];

    // Compute the LU factorization of a M by N matrix A
    cgetrf_(&n, &n, a, &lda, ipiv, &info);
    // Generate inverse of the matrix given its LU decompsotion
    cgetri_(&n, a, &lda, ipiv, works, &lwork, &info);
}

/*
 * Rotatating
 */

void transform(fcomplex *rho, fcomplex *x, fcomplex *rho_rotated, fcomplex *x_dag, int n)
{

     fcomplex work[N];
     fcomplex res[N*N];
     fcomplex alpha;
     fcomplex beta;
     beta.re = 0.0f;
     beta.im = 0.0f;
     alpha.re = 1.0f;
     alpha.im = 0.0f;
     int ipvi[N];
     // X^-1 A X = A (it overwrites it) 
     blas_matrix_multiplication(x_dag, rho, res, alpha, beta, n);
     print_matrix("px",n,n,res,n);
     print_matrix("inverse",n,n,x,n);
     blas_matrix_multiplication(x, res, rho, alpha, beta, n);
     return;
     
}

/*
 * Propagate
 */
void euler_propagate(fcomplex *rho, fcomplex *hamiltian, int n){
  
  return;
}

/*
 * Eigenvalue problem
 */

void diagonalization(fcomplex *a, fcomplex *w, fcomplex *vr, int n)
{
    int lda = N, ldvl = N, ldvr = N, info, lwork; /// Locals for cgeev

    // rwork dimension should be at least 2*n 
    float rwork[2*N];
    fcomplex wkopt;
    fcomplex* work;
    
    fcomplex vl[N*N];


    // Solve eigenproblem 
    lwork = -1;
    cgeev_( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
    &wkopt, &lwork, rwork, &info );
    printf("%f, %f , %d",wkopt.re,wkopt.im,info);
    lwork = (int)wkopt.re;
    //lwork = (int)6;
    work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
    cgeev_( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
    work, &lwork, rwork, &info );
    //print_matrix("hamiltonian",n,1,work,ldvl);
    if (info>0){
        printf("The algorithm failed to compute eigenvalues");
        exit(1);
    } 
    return;
}

/*
 * Record density of states
 */
void record(FILE *tape, fcomplex *rho, int step, int n){
     int i;
     fprintf(tape, "%d,", step);
     for(i=0; i< n; i++){
        fprintf(tape, "%f,", rho[IDX(i,i,n)].re);
     }
     fprintf(tape, "\n");
}

/*
 * Main
 */
int main( int argc, char *argv[] )
{
    /* Locals */
    int n = N;
    int i;
    fcomplex eigenvalues[N]; //eigenvalues values
    
    fcomplex eigenvectors[N*N];
    fcomplex eigenvectors_inv[N*N];

    fcomplex hamiltonian[6*6];
    fcomplex rho[N*N];
    fcomplex rho_rotated[N*N];
    /* Outputs */
    FILE *tape;
    tape = fopen("density.txt","a");


    /* get hamiltonian */
    get_hamiltonian(hamiltonian,n);
    diagonalization(hamiltonian, eigenvalues, eigenvectors, n);
    memcpy(eigenvectors_inv, eigenvectors,n);

    /* get X inv */
    inversematrix(eigenvectors_inv,n);

    /* initializing rho */
    get_initialrho(rho,n);
    record(tape,rho,0,n);

    /* rotating rho */
    transform(rho, eigenvectors, rho_rotated, eigenvectors_inv,n);
    print_matrix("transformed",n,n,rho,n);
 
    /* getting hamiltonian diag*/
    get_diaghamiltonian(hamiltonian,eigenvalues,n);

 
    /* propagation */
    for (i=0; i < STEPS; i++){
        propagate(rho_rotated,hamiltonian,n);
        transform(rho_rotated, eigenvectors_inv, rho, eigenvectors,n);
        record(tape,rho,i+1,n);
        
    } 
    
    return 0;
    
}


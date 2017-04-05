#ifndef _HYP2F1_H_
#define _HYP2F1_H_

#include <cmath>
#include <complex>
using namespace std;

void ggpgpp_exact( double t, double nu, double lambda, double Omega, double T, double kB, double hbar,
                   complex<double> &g, complex<double> &gp, complex<double> &gpp );

#define DBL_EPSILON       2.2204460492503131e-16    // float.h, 2^(-52)
#define DBL_MIN           2.2250738585072014e-308   // float.h, 2^(-1021)
//#define SQRT_DBL_EPSILON  1.4901161193847656e-08    // 2^(-26)
#define PI                3.14159265358979324
#define sqrt2Pi           2.50662827463100050       // = sqrt(2*Pi)

#define switchGAMMA       20  // switching between recursion and reflection

#define RADIUS            0.9 // 1843.0/2048.0 // ~ 0.8999... for hypgeom series
//#define STEPSIZE          3.0/1024.0    // ~ 0.00268... numerical differentiation

#define regRADIUS         1.0/1024.0    // for averaging around (numerical) poles

typedef std::complex<double> cplx;

#define ZERO              cplx(0.0, 0.0)
#define ONE               cplx(1.0, 0.0)
#define theI              cplx(0.0, 1.0)

// the 8 unit roots on the unitCircle * 4*regRADIUS, ignoring symmetries ...
static cplx ROOTS[] =
{
   cplx(+1.0                  * 4.0*regRADIUS, +0.0                 * 4.0*regRADIUS),
   cplx(+0.707106781186547525 * 4.0*regRADIUS, +.707106781186547525 * 4.0*regRADIUS),
   cplx(+0.0                  * 4.0*regRADIUS, +1.0                 * 4.0*regRADIUS),
   cplx(-0.707106781186547525 * 4.0*regRADIUS, +.707106781186547525 * 4.0*regRADIUS),
   cplx(-1.0                  * 4.0*regRADIUS, +0.0                 * 4.0*regRADIUS),
   cplx(-0.707106781186547525 * 4.0*regRADIUS, -.707106781186547525 * 4.0*regRADIUS),
   cplx(+0.0                  * 4.0*regRADIUS, -1.0                 * 4.0*regRADIUS),
   cplx(+0.707106781186547525 * 4.0*regRADIUS, -.707106781186547525 * 4.0*regRADIUS)
};
#define degROOTS  8 // = number of entries, should be chosen even

#define RE(z)     ((double)(std::real((cplx)z)))
#define IM(z)     ((double)(std::imag((cplx)z)))
#define ABS(z)    ((double)(std::abs((cplx)z)))
#define ARG(z)    ((double)(std::arg((cplx)z)))
//#define CONJUGATE(z) ((cplx)(std::conj((cplx)z)))

#define SIN(z)    ((cplx)(std::sin((cplx)z)))
#define SQRT(z)   ((cplx)(std::sqrt((cplx)z)))
#define EXP(z)    ((cplx)(std::exp((cplx)z)))
#define LOG(z)    ((cplx)(std::log((cplx)z)))
#define POW(z,a)  ((cplx)(std::pow((cplx)z,(cplx)a)))
//cplx POW(cplx z, cplx a); // this is NOT the one in the std library

// numTools
//int areClose (const double x1, const double x2, const double epsilon);
//int areEqual (cplx a, cplx b);

int roundToInteger (double x);
//int isNumericalInteger (cplx z);
int isNonNegativeInteger (cplx z);
int isOdd (int m);

double moduloPI(double x);
cplx rSIN (cplx mu);


// gamma
cplx pochhammer(cplx z, int n);
cplx GAMMA(cplx z);
cplx rGAMMA (cplx a);

// seriesCase
cplx series_2F1 (cplx A, cplx B, cplx C, cplx Z);
cplx polynomialseries_2F1 (cplx a, cplx b, cplx c, cplx z);
cplx polynom_2F1 (cplx a, cplx b, cplx c, cplx z,
                  int deg);

// Gosper
cplx Gosper_2F1 (cplx a, cplx b, cplx c, cplx z);

// linTransforms
cplx Euler (cplx a, cplx b, cplx c, cplx z);
cplx Euler_polynom (cplx a, cplx b, cplx c, cplx z);
cplx Pfaff (cplx a, cplx b, cplx c, cplx z);
cplx AS_15_3_6 (cplx a, cplx b, cplx c, cplx z);
cplx AS_15_3_7 (cplx a, cplx b, cplx c, cplx z);
//cplx AS_15_3_9_Leb (cplx a, cplx b, cplx c, cplx z);


// exceptionals
//cplx exceptional_15_3_6 (cplx a, cplx b, cplx c, cplx z);
//cplx exceptional_15_3_7 (cplx a, cplx b, cplx c, cplx z);
//cplx exceptional_15_3_9 (cplx a, cplx b, cplx c, cplx z);

// the actual hypergeometric functions 2F1
cplx reduce_2F1 (cplx a, cplx b, cplx c, cplx z);

cplx hyp2F1 (cplx a, cplx b, cplx c, cplx z);

// interface, structure for Excel
typedef struct
{
   double re,im;
}
cplx_structure;

#pragma fp_contract (off)

// round to nearest integer
int roundToInteger (const double x)
{
   if (0.0 < x)
      return (int)floor(x + 0.5);
   return (int)ceil(x - 0.5);
}

int isOdd (const int m)
{
   if ( (abs(m) % 2) == 0)
      return (0);
   else
      return (1);
}

// z = 0, 1, 2, ... ---> return 1, for me Naturals start in 0
int isNonNegativeInteger (const cplx z)
{
   double x = fabs(RE(z));

   if ( (IM(z) == 0.0) && (0.0 <= RE(z)) )
      if ( roundToInteger(x) == x )
         return(1);
   return(0);
}

// ############################################################################################

/*
 Cody & Waite
 returns x modulo Pi, ie: x = xi + k*PI, 0 <= |xi| < PI and k integer
 = redupi of Moshier's Cephes
 or GSL for Pi/2 (gsl_sf_sin_e) or also SUN's fdlibm
 There are certainly better range reductions (Payne&Hanek? Muller?)
 but it is exact except the last place for |x| <= 10^5, which is more
 than what is needed here and it is quite handy.
 */
double moduloPI(const double x)
{
   double DP1 = 3.14159265160560607910E0;
   double DP2 = 1.98418714791870343106E-9;
   double DP3 = 1.14423774522196636802E-17;
   double t;
   long k;

   if ( fabs(x) < PI) return x;

   t = x/PI;
   if( t >= 0.0 )
      t += 0.5;
   else
      t -= 0.5;

   k = (long)t;	/* the multiple */
   t = (double)k;
   t = ((x - t * DP1) - t * DP2) - t * DP3;
   // ensure that order of evaluation, the compiler should NOT optimize
   // that away to give t = x - k*(DP1+DP2+DP3) = x - k*Pi
   // fp_contract = off is used for that
   return (t);
}

// there *should* be no need to use this range reduction
cplx rSIN (const cplx z)  // = PI / sin( PI * z)
{
   cplx W_red=ZERO;
   double X=0.0, X_red=0.0;
   double x=0.0, y=0.0, t=0.0;
   int k=0;
   cplx result=ZERO;

   x = RE(z);
   y = IM(z);

   if ( x == 0.0 ) // sin(y*I) = sinh(y)*I, Pi/sin(Pi*y*I) = -I*Pi/sinh(Pi*y)
   {
      t = PI * y;
      result = - PI / sinh(PI * y);
      result = result * theI;
      return (result); // purely imaginary
   }

   // range reduction
   X = x*PI;
   X_red = moduloPI(X); // X = X_red + k*Pi
   k = roundToInteger( (X - X_red)/PI );

   W_red = X_red + y*PI *theI;

   result = PI / SIN(W_red);
   if ( isOdd(k) == 1 ) //if ( isEven(m) == 1 )
      result = - result;

   if (y == 0.0 )  // purely real
      result = RE(result) + 0.0*theI;
   return (result);
}

/*
 This computes Gamma(z) and the reciprocal 1/Gamma(z), z complex.

 Essentially it is just Godfrey's implementation of the Lanczos method.
 However I (besides some pre-computed values) added some tweaks around
 (negative) integers, so it is longer and looks more complicated. The
 reason for me was, that in those points (real part < 0) I got results,
 which seemed not precise enough for my needs. However I am not sure,
 whether that was related to the MS compiler or my coding style (but
 looking in GSL sources indicates others had that problem as well).
 */

// ********************************************************************
// Gamma(z)
// ********************************************************************

// pochhammer(z,n) = `z*(z+1)*...*(z+n-1)` = GAMMA(n+z)/GAMMA(z)
// this is a brute version
cplx pochhammer(const cplx z, const int n)
{
   int k=0;
   cplx p;

   if ( z == ZERO && n == 0 )
      return ONE;
   if ( z == ZERO )  // n <> 0
      return ZERO;
   if ( n < 0 )      // write as Gamma, to see it is true
      return ONE/pochhammer(z + (cplx)n, -n);
   if ( n == 0 )
      return ONE;
   if ( n == 1 )
      return z;
   if ( n == 2 )
      return z * (z + ONE);

   p = z;
   for(k=1; k <= n-1; k++)
      p = p * (z + (cplx)k);

   return (cplx)p;
}

// = (1/(z+1) + GAMMA(z)) - 1/(z+1) in z = -1
// series around -1 up to | z - (-1) | <= 1/128
cplx G_minus1(const cplx z)
{
   cplx w = z + ONE;   // close to 0
   cplx w2 = w * w;
   cplx w4 = w2 * w2;
   cplx s = ZERO;

   // sum up starting with the smallest term
   // do not use a Horner form in 0 because of very small w
   s =  -0.501241626455653050 * w2*w4;
   s += -1.49724338689808458  * w*w4;
   s += -0.504094272276808391 * w4;
   s += -1.48608934117995359  * w*w2;
   s += -0.504361254345553406 * w2;
   s += -1.41184033042643969  * w;
   s += -0.422784335098467139;

   s = s - ONE/w;
   return s;
}

// = (Gamma(z) - 1/z) + 1/z
cplx G_minus0(const cplx z)
{
   cplx z2 = z * z;
   cplx z4 = z2 * z2;
   cplx s = ZERO;

   s =  -0.996001760442431534  * z2*z4 ;
   s += +0.993149114621276193  * z*z4;
   s += -0.981995068903145202  * z4;
   s += +0.981728086834400187  * z*z2;
   s += -0.907479076080886289  * z2;
   s += +0.989055995327972555  * z;
   s += -0.577215664901532861;

   s = s + ONE/z;
   return s;
}

// following Godfrey, see his Matlab implementation gamma.m
cplx GAMMA(const cplx Z)
{
   cplx z;
   double x;
   cplx zh, zgh, zp;
   cplx sum;
   int k=0;
   cplx result;
   int nx = 0; // integer part of real(z)
   double xi=0.0;

   static double g = 607.0/128.0;
   static double c[] =
   {
      0.99999999999999709182,
      57.156235665862923517,
      -59.597960355475491248,
      14.136097974741747174,
      -0.49191381609762019978,
      0.33994649984811888699e-4,
      0.46523628927048575665e-4,
      -0.98374475304879564677e-4,
      0.15808870322491248884e-3,
      -0.21026444172410488319e-3,
      0.21743961811521264320e-3,
      -0.16431810653676389022e-3,
      0.84418223983852743293e-4,
      -0.26190838401581408670e-4,
      0.36899182659531622704e-5
   };
   int dimc = 15; // length of c

   // exact values for GAMMA(n), n = 0 ... 19
   static double G[] =
   {
      1.0000000000000000e+308, // to indicate overflow
      1.0,
      1.0,
      2.0,
      6.0,
      24.0,
      120.0,
      720.0,
      5040.0,
      40320.0,
      362880.0,
      3628800.0,
      39916800.0,
      479001600.0,
      6.2270208000000000e+09,
      8.7178291200000000e+10,
      1.3076743680000000e+12,
      2.0922789888000000e+13,
      3.5568742809600000e+14,
      6.4023737057280000e+15,
      1.2164510040883200e+17
   };

   z = (cplx)Z;
   x = RE(z);
   nx = roundToInteger( RE(z) );
   xi = x - nx;


   if ( ( IM(z)== 0.0) && (x == int(x)) )
   {
      if (nx <= 0) return (cplx(1e308, 0.0));   // indicating pole for z =
      if (nx <= 19) return (cplx)G[nx];         // return exact value, if tabled
   }

   if ( nx < -switchGAMMA )
      return rSIN(z) / GAMMA(-z + ONE);         // reflection for small inputs

   if ( nx == 0 )
   {
      result = (cplx)ONE/pochhammer(z, 1);
      result = result * GAMMA(z + ONE);
      return result;
   }

   if ( nx == -1 && ABS(z + ONE) <= 1.0/128.0 ) // small circle around -1
   {
      result = G_minus1(z);
      return result;
   }

   /* not needed
    if ( nx == -2 && ABS(z + ONE+ONE) <= 2.0/128.0 ) // small circle around -2
    {
    result = G_minus2(z);
    return result;
    }
    */
   if ( nx <= -1 && ABS(z - (cplx)nx ) <= 1.0/128.0 )
   {
      result = (cplx)ONE/pochhammer(z, -nx - 1);
      result = result * GAMMA(z - (cplx)nx - ONE);
      return result;
   }

   if ( nx <= -1 )
   {
      result = (cplx)ONE/pochhammer(z, -nx + 1);
      result = result * GAMMA(z - (cplx)nx + ONE);
      return result;
   }

   // Lanczos method for Gamma, according to Godfrey
   z = z - ONE;

   // his trick for avoiding FP overflow above z=141
   zh = z + 0.5*ONE;
   zgh=zh + (cplx)g;
   zp=POW(zgh, (zh*0.5));

   // starting with maximal index and descending, so with increasing summands
   sum = ZERO;
   for (k=dimc - 1; 0<k; k--)
   {
      sum = sum + (cplx)c[k]/(z + (cplx)k);
   }

   result = sqrt2Pi * (sum + (cplx)c[0]); // final, minimal index
   result = result * zp*zp * EXP(-zgh);
   return result;
}

// ********************************************************************
// 1 / Gamma(z)
// ********************************************************************

// series for 1/Gamma(z) around z=0, radius 1/128
cplx rG_minus0 (const cplx z)
{
   cplx result;
   cplx z2 = z * z;
   cplx z3 = z * z2;
   cplx z4 = z2 * z2;
   cplx s = ZERO;

   s =  +0.721894324666309954e-2 * z4*z4;
   s += -0.962197152787697356e-2 * z3*z4;
   s += -0.421977345555443367e-1 * z3*z3;
   s += +0.166538611382291490    * z2*z3;
   s += -0.420026350340952355e-1 * z4;
   s += -0.655878071520253881    * z3;
   s += +0.577215664901532861    * z2;
   s += z;

   return (s);
}

// series for 1/Gamma(z) around z=-1, radius 1/128
cplx rG_minus1 (cplx z)
{
   cplx w = z + ONE;   // close to 0
   cplx w2 = w * w;
   cplx w3 = w * w2;
   cplx w4 = w2 * w2;
   cplx s = ZERO;

   s =  -0.168409147745400731e-1 * w4*w4;
   s += -0.325757630276673632e-1 * w3*w4;
   s += +0.208736345937835826    * w3*w3;
   s += -0.208541246416386725    * w2*w3;
   s += -0.613875436486158646    * w4;
   s += +1.23309373642178674     * w3;
   s += +0.422784335098467139    * w2;
   s += - w;

   return s;
}

// following Godfrey, see his Matlab implementation gamma.m
cplx rGAMMA(const cplx Z)
{
   cplx z;
   cplx z_red;
   double x;
   cplx zh, zgh, zp;
   cplx rzp;
   cplx sum;
   int k=0;
   cplx result;
   int nx = 0; // integer part of real(z)
   cplx W_red=ZERO;
   double X=0.0, X_red=0.0;

   static double g = 607.0/128.0;
   static double c[] =
   {
      0.99999999999999709182,
      57.156235665862923517,
      -59.597960355475491248,
      14.136097974741747174,
      -0.49191381609762019978,
      0.33994649984811888699e-4,
      0.46523628927048575665e-4,
      -0.98374475304879564677e-4,
      0.15808870322491248884e-3,
      -0.21026444172410488319e-3,
      0.21743961811521264320e-3,
      -0.16431810653676389022e-3,
      0.84418223983852743293e-4,
      -0.26190838401581408670e-4,
      0.36899182659531622704e-5
   };
   int dimc = 15; // length of c

   // exact values in (positive) integers, 1/Gamma(n), n = 0 ... 19
   static double rG[] =
   {
      0.0,
      1.0000000000000000e+00,
      1.0000000000000000e+00,
      5.0000000000000000e-01,
      1.6666666666666667e-01,
      4.1666666666666667e-02,
      8.3333333333333333e-03,
      1.3888888888888889e-03,
      1.9841269841269841e-04,
      2.4801587301587302e-05,
      2.7557319223985891e-06,
      2.7557319223985891e-07,
      2.5052108385441719e-08,
      2.0876756987868099e-09,
      1.6059043836821615e-10,
      1.1470745597729725e-11,
      7.6471637318198165e-13,
      4.7794773323873853e-14,
      2.8114572543455208e-15,
      1.5619206968586226e-16,
      8.2206352466243297e-18
   };

   z = (cplx)Z;
   x = RE(z);
   nx = roundToInteger( RE(z) );

   if ( ( IM(z)== 0.0) && (x == int(x)) )
   {
      if (nx <= 0) return ZERO;
      if (nx <= 19) return (cplx)rG[nx];      // return exact value, if tabled
   }

   if( nx < -switchGAMMA ) // reflection for small inputs, argument reduction
   {
      X = x*PI;
      X_red = moduloPI(X); // X = X_red + k*Pi
      nx = roundToInteger( (X - X_red)/PI );
      W_red = X_red + IM(z)*PI *theI;
      result = PI / SIN(W_red)/PI;
      if (IM(z) == 0.0 )
         result = RE(result);
      if ( isOdd(nx) == 1 )
         result = - result;      // sin(W_red +-n*Pi) = sin(W_red)*(-1)^n
      result = (cplx)result * (cplx)GAMMA(ONE - z);
      return result;
   }

   if ( nx == 0 && ABS(z) <= 1.0/128.0 ) // small circle around 0
      return rG_minus0(z);

   if ( nx == -1 && ABS(z + ONE) <= 1.0/128.0 ) // small circle around -1
      return rG_minus1(z);

   // not needed ?
   //if ( nx == -2 && ABS(z + ONE+ONE) <= 2.0/128.0 ) // small circle around -2
   //return rG_minus2(z);

   // 1/GAMMA(z) = pochhammer(z,k)/GAMMA(z+k), k = -nx - 1, so Re(z+k) ~ -1
   if ( nx <= -1 && ABS(z - (cplx)nx ) <= 1.0/128.0 )
   {
      result = (cplx)pochhammer(z, -nx - 1);
      result = result / GAMMA(z - (cplx)nx - ONE);
      return result;
   }

   // 1/GAMMA(z) = pochhammer(z,k)/GAMMA(z+k), k = -nx + 1, so Re(z+k) ~ 1
   if ( nx <= -1 )
   {
      result = (cplx)pochhammer(z, -nx + 1);
      result = result / GAMMA(z - (cplx)nx + ONE);
      return result;
   }

   // Lanczos method for Gamma, according to Godfrey
   z = z - ONE;

   // his trick for avoiding FP overflow above z=141
   zh = z + 0.5*ONE;
   zgh=zh + (cplx)g;
   zp=POW(zgh, (zh*0.5));

   // starting with maximal index and descending, so with increasing summands
   sum = ZERO;
   for (k=dimc - 1; 0<k; k--)
   {
      sum = sum + (cplx)c[k]/(z + (cplx)k);
   }

   result = sqrt2Pi * (sum + (cplx)c[0]); // final, minimal index
   //result = result * zp*zp * EXP(-zgh);
   rzp = POW(zgh, -(zh*0.5));
   result = ONE/result * rzp*rzp * EXP(zgh);
   return result;
}

// Euler's transform, no restrictions for parameters or variable
cplx Euler (const cplx a, const cplx b, const cplx c, cplx z)     // = A&S 15.3.3
{
   return POW(ONE - z, c-a-b) * (cplx)reduce_2F1(c - b, c - a, c, z);
}

// polynomial case only to be used, if c-b or c-a are integers <= 0
cplx Euler_polynom (const cplx a, const cplx b, const cplx c, cplx z)
{
   return POW(ONE - z, c - a - b) *
          (cplx)polynomialseries_2F1(c - b, c - a, c, z);
}

// Pfaff's transform, (1-z)^(-a)*hypergeom([a, c-b],[c],z/(z-1)),  w = z/(z-1)
// no restrictions for parameters, but z off from 1, A&S 15.3.4 or 15.3.5
cplx Pfaff (const cplx a, const cplx b, const cplx c, const cplx z)
{
   return POW(ONE - z, - a) * (cplx)reduce_2F1(a, c - b, c, z/(z-ONE));
}

//############################################################################

/*
 The main linear transforms from Abramowitz & Stegun.

 Note that this can lead to so called  "catastrophic cancellations" errors
 (same magnitude, opposite sign, hence only few final decimal places survive,
 which is not enough for precise results or even quite false in sequel).
 That is independent from prior exact results and just depending on input.

 The most reasonably way is to use a few as possible transform, having just
 a single one is desirably.

 Also observe, that for c ~ negative integer this may become almost useless
 by the poles of GAMMA, especially if c becomes smaller.

 For small 'exceptional' cases Cauchy's formula is used.
 */


cplx AS_15_3_6 (const cplx a, const cplx b, const cplx c, cplx z)   // w = 1 - z
{
   cplx w = (ONE - z);
   cplx mu = c - a - b;
   cplx f1=ZERO, f2=ZERO, g1=ZERO, g2=ZERO;
   double m = (double)roundToInteger(RE(mu));
   cplx center=a;
   cplx s=ZERO;
   int k=0;

   // close to exceptional: Cauchy's formula, discrete version, the
   // function is holomorphic in the 1st parameter (if others are fixed)
   if ( ABS(mu - cplx(m,0.0) ) < regRADIUS ) // if ( ABS(c-a-b - m ) < regRADIUS )
   {
      for(k=0; k<degROOTS/2; k++) // the other roots are the negative ones
         s = s +
             AS_15_3_6(+ROOTS[k] + center,b,c, z) +
             AS_15_3_6(-ROOTS[k] + center,b,c, z);

      s = 1.0/(double)degROOTS *s; // average around center
      return s;
   }

   f1 = GAMMA(mu) * rGAMMA(b + mu) * rGAMMA(a + mu);
   f2 = POW(w, mu) * GAMMA(-mu) * rGAMMA(a) * rGAMMA(b);

   g1 = reduce_2F1(a, b, -mu + ONE, w);
   g2 = reduce_2F1(a + mu, b + mu, ONE + mu, w);
   return GAMMA(c)*(f1 * g1 + f2 * g2);
}

// almost identical, just other functions as input
cplx AS_15_3_7 (const cplx a, const cplx b, const cplx c, cplx z)   // w = 1/z
{
   cplx w = (ONE / z);
   cplx mu = b-a;
   cplx f1=ZERO, f2=ZERO, g1=ZERO, g2=ZERO;
   double m = (double)roundToInteger(RE(mu));
   cplx center=a;
   cplx s=ZERO;
   int k=0;

   if ( ABS(mu - cplx(m, 0.0) ) < regRADIUS )
   {
      for(k=0; k<degROOTS/2; k++)
         s = s +
             AS_15_3_7(+ROOTS[k] + center,b,c, z) +
             AS_15_3_7(-ROOTS[k] + center,b,c, z);

      s = 1.0/(double)degROOTS *s;
      return s;
   }

   f1 = POW(-z, -a) * GAMMA(mu) * rGAMMA(c - a) * rGAMMA(b);
   f2 = POW(-z, -b) * GAMMA(-mu) * rGAMMA(a) * rGAMMA(c - b);

   g1 = reduce_2F1(a, ONE + a - c, -mu + ONE, w);
   g2 = reduce_2F1(b, b - c + ONE, mu + ONE, w);
   return GAMMA(c) * ( f1 * g1 + f2 *g2 );
}

cplx branchPoint_2F1 (cplx a, cplx b, cplx c)
{
   // branching point z = 1:
   //if (0.0 == IM(z) && 1.0 == RE(z))
   {
      //care only for the classical case
      if ( 0.0 < RE(c - a - b) )
         return(GAMMA(c) * GAMMA(c - a - b) * rGAMMA(c - a) * rGAMMA(c - b));
      // otherwise imitate limit from the left by smallest IEEE number < 1
      hyp2F1(a,b,c, cplx(1.0 - DBL_EPSILON*0.5, 0.0));
   }
   return ZERO;
}

cplx hyp2F1 (const cplx a, const cplx b, const cplx c, const cplx Z)
{
   double x=RE(Z),y=IM(Z);
   cplx z=ZERO;
   cplx result;

   // branch cut
   if ( (0.0 == y) && (1.0 < x) )
      y = - (1e-306);
   z = x + y*theI;

   result = reduce_2F1(a,b,c,z);

   // clean up branch cut
   if ( fabs(IM(result)) <= 1e-306)
      result = RE(result) + ZERO;

   // the disk, real parameters give real results
   if ( (x < 1.0) && (0.0 == y)
         &&(IM(a) == 0.0) && (IM(b) == 0.0) && (IM(c) == 0.0) )
      result = RE(result) + ZERO;
   return result;
}

cplx reduce_2F1 (const cplx a, const cplx b, const cplx c, cplx z)
{
   cplx result;
   double delta;
   cplx w;
   cplx M;

   delta = 0.51776695296636881e-1; // = (sqrt(2.0) - 1.0) / 8.0

   if ( ZERO == z)
      return(ONE);

   if (isNonNegativeInteger(-c) == 1)
      return(0.1e308);

   if ( ABS(b) < ABS(a) )
      return reduce_2F1(b,a,c,z);

   // try to find good constellations for parameters by Euler's transform
   // note: applying twice gives identity (---> infinite loop by recursion),
   // so care for strict conditions to apply only once
   if ( ABS((c-a)*(c-b)/c) < ABS(a*b/c) )
      return Euler(a,b,c, z); // w = z, it is not changed

   if (isNonNegativeInteger(-a) == 1 || isNonNegativeInteger(-b) == 1)
      return(polynomialseries_2F1(a, b, c, z));

   if (isNonNegativeInteger(-(c - a)) == 1 || isNonNegativeInteger(-(c - b)) == 1)
      return(Euler_polynom(a, b, c, z)); // Euler's relation gives a polynom for that

   // reduction to some other elementary function, deactivate it for tests,
   // since test function often rely on that relations
   // power: (1-x)^(-a) = hypergeom([a],[],x)
   // log: ln(1-x) = -x*hypergeom([1, 1],[2],x)
   // asin: arcsin(x^(1/2)) = x^(1/2)*hypergeom([1/2, 1/2],[3/2],x)
   // atan: arctan(x^(1/2)) = x^(1/2)*hypergeom([1/2, 1],[3/2],-x)
   /*
    if (b == c) return POW(ONE - z, -a);
    if (a == c) return POW(ONE - z, -b);
    if (ONE == a && ONE == b && c == ONE+ONE) return -log(ONE - z)/z;
    */

   // branching point z = 1:
   if (0.0 == IM(z) && 1.0 == RE(z))
      return  branchPoint_2F1(a,b,c);

   // branch cut: from below
   if ( 0.0 == IM(z) && 1.0 < RE(z) )
      z = RE(z) - 1e-306 * theI;

   // special case: z ~ 1
   if ( ABS(z - ONE) <= 1.0/1024.0 && (ZERO != z) )
      return(AS_15_3_6(a, b, c, z));  // w = 1 - z ~ 0

   if ( ABS(z) <= RADIUS )
      return series_2F1(a,b,c,z);

   //if ( ABS(z - ONE) <= 1.0/1024.0 && (ZERO != z) )
   //return(AS_15_3_6(a, b, c, z));  // w = 1 - z;

   // = AS_15_3_8 = AS_15_3_7 @ Pfaff, w = 1 / [z/(z-1)] = 1/(z-1)
   if ( 0.0 <= RE(z) &&
         RADIUS < ABS(z) && ABS(z) <= 2.0 &&
         fabs(ARG(z)) <= 1.0/3.0 ) // && 1.0/RADIUS < ABS(z/(z-ONE)) )
      if ( ABS(a*(c - b)) < ABS((c - a)*b) )            // use A&S 15.3.4
         return POW(ONE - z, - a) * AS_15_3_7(a, c - b, c, z/(z-ONE));
      else                                              // use A&S 15.3.5
         return POW(ONE - z, - b) * AS_15_3_7(c - a, b, c, z/(z-ONE));

   if ( 2.0 <= ABS(z) )
      return AS_15_3_7(a,b,c,z);      // w = 1 / z;

   return Gosper_2F1(a, b, c, z);
}

cplx polynom_2F1 (cplx a, cplx b, cplx c, cplx z,
                  int deg)
{
   int n;
   int nMax;
   cplx r;
   cplx s;
   cplx t;

   if (2000 < deg)
      nMax = 2000;
   else
      nMax = deg;

   if (nMax == 0)
      return(ONE);
   if (ZERO == a || ZERO == b)
      return(ONE);
   if (ZERO == z)
      return(ONE);

   if (nMax == 1)
      return(ONE + a * b / c * z);

   n = 0;
   s = ONE;
   t = ONE;
   while (n < nMax)
   {
      n = n + 1;
      r = a * b / c / (cplx) n * z;
      if (ZERO == r)
         break;

      t = t * r;
      s = t + s;
      a = a + ONE;
      b = b + ONE;
      c = c + ONE;
   }

   // clean up for the real case
   if ( (IM(a) == 0.0) && (IM(b) == 0.0) && (IM(c) == 0.0) &&
         (IM(z) == 0.0) )
      s = RE(s) + 0.0*theI;

   return(s);
}

// to be used, if a or b are negative integers, then it is a polynom
// and one has to determine the degree
cplx polynomialseries_2F1 (const cplx a, const cplx b, const cplx c, const cplx z)
{
   int n=0;
   int int_a=0, int_b=0;
   cplx result=ZERO;
   double re_a = RE(a), re_b = RE(b);

   int_a = roundToInteger( re_a );
   int_b = roundToInteger( re_b );

   // a integer and <= 0
   if ( int_a <= 0 && IM(a) == 0.0 && re_a == (double)int_a )
      n = int_a;    // i.e. a = n
   // b integer and <= 0
   if ( int_b <= 0 && IM(b) == 0.0 && re_b == (double)int_b
         && int_b <= n )
      n = int_b;    // b = n
   // the smaller one is taken, n <= 0

   result = polynom_2F1(a, b, c, z, -n);
   return(result);
}

// hypergeometric series for abs(z) <= RADIUS < 1
cplx series_2F1 (cplx a, cplx b, cplx c, cplx z)
{
   int n;
   int nMax;
   cplx r;
   cplx s;
   cplx sNew;
   cplx t;
   int lastEqual = 0;

   nMax = 2000; // this is never reached

   if ( a == ZERO || b == ZERO )
      return(ONE);
   if (z == ZERO)
      return(ONE);

   n = 0;
   s = ONE;
   t = ONE;
   while (n < nMax)
   {
      n = n + 1;
      r = ( a * (b/c) ) / (cplx) n * z;
      t = (t * r);
      sNew = (t + s);

      if( s == sNew )
         if ( 1 == lastEqual ) break;
         else lastEqual = 1;
      else
         lastEqual = 0;

      s = sNew;

      a = (a + ONE);
      b = (b + ONE);
      c = (c + ONE);
   }

   // clean up
   if ( (IM(a) == 0.0) && (IM(b) == 0.0) && (IM(c) == 0.0) &&
         (IM(z) == 0.0) )
      s = RE(s) + 0.0*theI;

   return(s);
}

// Gosper's method
// http://www.math.utexas.edu/pipermail/maxima/2006/000126.html
// essential recursion is giveb below
cplx Gosper_2F1 (cplx a, cplx b, cplx c, cplx z)
{
   int k;
   cplx K;
   int kMax;
   cplx d_k;
   cplx e_k;
   cplx f_k;
   cplx d_k1;
   cplx e_k1;
   cplx f_k1, t;
   cplx xi;
   cplx lambda;
   cplx mu;
   cplx result;
   int lastEqual = 0;

   if (ABS(z) <= RADIUS)
   {
      result = series_2F1(a, b, c, z);
      return(result);
   }

   kMax = 2000; // never reached
   xi = z / (z - ONE);
   mu = c - a - b;
   d_k = ZERO;
   e_k = ONE;
   f_k = ZERO;
   for (k = 0; k <= kMax; k++)
   {
      K = cplx(double(k), 0.0);
      lambda = (K + a) * (K + b) / (K + ONE) / (2.0 * K + c) / (2.0 * K + c + ONE);
      d_k1 = (lambda * z * (xi * d_k * (mu + K) + e_k));
      e_k1 = (lambda * z * (-a * b * xi * d_k + (K + c) * e_k));
      t    =  d_k * (-K / (-ONE + z) - (b * a - (mu + K) * K) / (2.0 * K + c) * xi) + e_k;
      f_k1 = (f_k + t);

      if( f_k1 == f_k )
         if ( 1 == lastEqual ) break;
         else lastEqual = 1;
      else
         lastEqual = 0;

      d_k = (d_k1);
      e_k = (e_k1);
      f_k = (f_k1);
   }
   result = f_k;

   if ( (ABS(z) < 1.0) &&
         (IM(a) == 0.0) && (IM(b) == 0.0) && (IM(c) == 0.0) && (IM(z) == 0.0) )
      result = RE(result) + 0.0*theI;

   return(result);
}

static complex<double> digamma( complex<double> z )
{
   gsl_sf_result glnr,garg;
   // FIXME check definition!
   gsl_sf_complex_psi_e(real(z),imag(z),&glnr,&garg);
   // printf("%f %f\n",glnr.val,garg.val);
   return complex<double>(glnr.val,garg.val);
}

static complex<double> lerchphi( complex<double> z, complex<double> b, complex<double> v )
{
   //   fprintf(stderr,"lerchphi(%f %f|%f %f|%f %f)\n",real(z),imag(z),real(b),imag(b),real(v),imag(v));
   if (b!=1.0)
   {
      fprintf(stderr,"lerchphi domain error\n");
      exit(1);
   }
   return hyp2F1(1.0,v,1.0+v,z)/v;
}

static complex<double> HarmonicNumber( complex<double> z )
{
   return 0.57721566490153286061+digamma(1.0+z);
}

static double Pi=M_PI;

static complex<double> cot( complex<double> x )
{
   return 1.0/tan(x);
}

void ggpgpp_exact( double t, double nu, double lambda, double Omega, double T, double kB, double hbar,
                   complex<double> &g, complex<double> &gp, complex<double> &gpp )
{
   if(t==0)
   {
      g=0.0;
      gp=0.0;
      gpp=0.0;
   }
   else
   {
      complex<double> BB,CC,DD,EE,FF,GG,II,JJ,KK,NN,XX,YY,ZZ;

      BB=HarmonicNumber(-(hbar*nu)/(2.*kB*Pi*T) - ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      CC=HarmonicNumber(+(hbar*nu)/(2.*kB*Pi*T) - ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      DD=HarmonicNumber(-(hbar*nu)/(2.*kB*Pi*T) + ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      EE=HarmonicNumber(+(hbar*nu)/(2.*kB*Pi*T) + ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));

      JJ=cot((hbar*nu)/(2.*kB*T) - ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*T));
      KK=cot((hbar*nu)/(2.*kB*T) + ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*T));

      XX=exp(-nu*t - (complex<double>(0.0,1.0))*Omega*t);
      YY=exp(-nu*t + (complex<double>(0.0,1.0))*Omega*t);
      ZZ=expl((-2*kB*Pi*t*T)/hbar);

      FF=lerchphi(ZZ,1,1 - (hbar*nu)/(2.*kB*Pi*T) - ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      GG=lerchphi(ZZ,1,1 + (hbar*nu)/(2.*kB*Pi*T) - ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      II=lerchphi(ZZ,1,1 - (hbar*nu)/(2.*kB*Pi*T) + ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));
      NN=lerchphi(ZZ,1,1 + (hbar*nu)/(2.*kB*Pi*T) + ((complex<double>(0.0,0.5))*hbar*Omega)/(kB*Pi*T));

      g=
         (lambda*(DD*nu + EE*nu + cplx(0,1)*DD*Omega -
                  cplx(0,1)*EE*Omega + cplx(0,2)*nu*Pi -
                  JJ*nu*Pi - KK*nu*Pi - cplx(0,1)*JJ*Omega*Pi +
                  cplx(0,1)*KK*Omega*Pi - DD*pow(nu,2)*t +
                  EE*pow(nu,2)*t - DD*pow(Omega,2)*t +
                  EE*pow(Omega,2)*t - cplx(0,2)*pow(nu,2)*Pi*t +
                  JJ*pow(nu,2)*Pi*t + KK*pow(nu,2)*Pi*t -
                  cplx(0,2)*pow(Omega,2)*Pi*t +
                  JJ*pow(Omega,2)*Pi*t + KK*pow(Omega,2)*Pi*t +
                  BB*(nu - cplx(0,1)*Omega -
                      (pow(nu,2) + pow(Omega,2))*t) +
                  CC*(nu + pow(nu,2)*t +
                      Omega*(cplx(0,1) + Omega*t)) -
                  cplx(0,1)*nu*Pi*XX + KK*nu*Pi*XX - Omega*Pi*XX -
                  cplx(0,1)*KK*Omega*Pi*XX - cplx(0,1)*nu*Pi*YY +
                  JJ*nu*Pi*YY + Omega*Pi*YY +
                  cplx(0,1)*JJ*Omega*Pi*YY +
                  ((FF + GG + II + NN)*nu -
                   cplx(0,1)*(FF - GG - II + NN)*Omega)*ZZ +
                  4*nu*log(1.0 - ZZ)))/
         (2.*(pow(nu,2) + pow(Omega,2))*Pi);

      gp=(lambda*(-BB + CC - DD + EE + Pi*(JJ + KK - KK*XX - JJ*YY
                                           + cplx(0,1)*(-2.0 + XX + YY)) + (-FF + GG - II + NN)*ZZ))/(2.*Pi);

      gpp=(lambda*Pi*((cplx(0,-1) + KK)*(nu + cplx(0,1)*Omega)*XX +
                      (cplx(0,-1) + JJ)*(nu - cplx(0,1)*Omega)*YY) +
           lambda*((FF + GG + II + NN)*nu + cplx(0,1)*(FF - GG - II + NN)*Omega)*ZZ)/(2.*Pi);
   }
}

static void ggpgpp_for_sd( double t, double T, int Npeaks, double* lambda, double_complex* gamma,
                           double_complex &g, double_complex &gp, double_complex &gpp)
{
   double_complex rg,rgp,rgpp;

   int k=0;
   if(Npeaks==1 and imag(gamma[k])==0.0)
   {
      ggpgpp_exact(t,real(gamma[k]),lambda[k],imag(gamma[k]),T,const_kb,const_hbar,rg,rgp,rgpp);
      g+=rg;
      gp+=rgp;
      gpp+=rgpp;
   }
   else
   {
      for(int k=0; k<Npeaks; k++)
      {
         ggpgpp_exact(t,real(gamma[2*k+0]),lambda[2*k+0],imag(gamma[2*k+0]),T,const_kb,const_hbar,rg,rgp,rgpp);
         g+=rg;
         gp+=rgp;
         gpp+=rgpp;
         ggpgpp_exact(t,real(gamma[2*k+1]),lambda[2*k+1],imag(gamma[2*k+1]),T,const_kb,const_hbar,rg,rgp,rgpp);
         g+=rg;
         gp+=rgp;
         gpp+=rgpp;
      }
   }
}

// this is required to pass parameters to static non-class member
// function for GSL integrator
struct f_params
{
   // must pass everything
   bool modified_redfield; // true: return modified RF, false: return generlized Foerster
   int ek;
   int ekp;
   long double Ek;
   long double Ekp;
   long double *Evk;
   long double *Evkp;
   double* l_list;
   double Temp;
   double** sd_lambda_invcm;
   double** sd_gamma_fs;
   double** sd_gammapos_invcm;
   vector<int> sd_CoupledSite;
   vector<int> sd_Npeaks;
   // needed for combined method
   long double HCW;
   int part;

   // write some results back to this variable
   double re_gpkkkkp;
   double lkkkk;
};

typedef long double long_double;
typedef complex<long double> long_double_complex;

double f(double t, void * p)
{
   f_params &params= *reinterpret_cast<f_params *>(p);

   double *lambda;
   double_complex* gamma;

   long_double lkkkk=0.0;
   long_double lkpkpkpkp=0.0;
   long_double lkpkpkk=0.0;
   long_double lkpkkpkp=0.0;
   long_double lkpkpkkp=0.0;
   long_double_complex gkkkk=0.0;
   long_double_complex gkpkpkpkp=0.0;
   long_double_complex gkpkpkk=0.0;
   long_double_complex gpkpkkpkp=0.0;
   long_double_complex gpkpkkk=0.0;
   long_double_complex gpkpkpkkp=0.0;
   long_double_complex gpkkkkp=0.0;
   long_double_complex gppkkpkpk=0.0;

   // sum over all coupled sites
   for(int i=0; i<params.sd_CoupledSite.size(); i++)
   {
      // obtain lambda and gamma for this site
      if(params.sd_Npeaks[i]==1 and params.sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[params.sd_Npeaks[i]];
         gamma=new double_complex[params.sd_Npeaks[i]];
         for(int k=0; k<params.sd_Npeaks[i]; k++)
         {
            lambda[k]=params.sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
      }
      else
      {
         lambda=new double[2*params.sd_Npeaks[i]];
         gamma=new double_complex[2*params.sd_Npeaks[i]];
         for(int k=0; k<params.sd_Npeaks[i]; k++)
         {
            lambda[2*k]=params.sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            lambda[2*k+1]=params.sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            gamma[2*k]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15),+params.sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
            gamma[2*k+1]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15),-params.sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
      }

      // compute g,g',g'' at this site at time t (needs to know lambda, gamma, Omega)

      double_complex g,gp,gpp;
      long_double_complex ldg,ldgp,ldgpp;

      ggpgpp_for_sd(t,params.Temp,params.sd_Npeaks[i],lambda,gamma,g,gp,gpp);
      // fprintf(stderr,"t=%e g=(%e,) gp=(%e,) gpp=(%e,)\n",t,real(g),real(gp),real(gpp));
      ldg=g;
      ldgp=gp;
      ldgpp=gpp;

      delete[] lambda;
      delete[] gamma;

      // compute g,g',g'' at this site at time t (needs to know lambda, gamma, Omega)

      int siteindex=params.sd_CoupledSite[i];
      gkkkk    +=((long_double)powl(params.Evk [siteindex],4)*ldg);
      gkpkpkpkp+=((long_double)powl(params.Evkp[siteindex],4)*ldg);
      gkpkpkk  +=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*ldg;

      gpkpkkpkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;
      gpkpkkk  +=long_double(params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*ldgp;
      gpkpkpkkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;
      gpkkkkp  +=long_double(params.Evk [siteindex]*params.Evk [siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;

      gppkkpkpk+=long_double(params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex])*ldgpp;

      lkkkk    +=((long_double)powl(params.Evk [siteindex],4)*params.l_list[siteindex]);
      lkpkpkpkp+=((long_double)powl(params.Evkp[siteindex],4)*params.l_list[siteindex]);

      lkpkpkk +=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*params.l_list[siteindex];
      lkpkkpkp+=long_double(params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex])*params.l_list[siteindex];
      lkpkpkkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*params.l_list[siteindex];
   }
   gkkkk/=const_hbar;

   gkpkpkpkp/=const_hbar;
   gkpkpkk/=const_hbar;
   gpkpkkpkp/=const_hbar;
   gpkpkkk/=const_hbar;
   gpkpkpkkp/=const_hbar;
   gpkkkkp/=const_hbar;

   gppkkpkpk/=const_hbar;

   lkkkk/=const_hbar;
   lkpkpkpkp/=const_hbar;
   lkpkpkk/=const_hbar;
   lkpkkpkp/=const_hbar;
   lkpkpkkp/=const_hbar;

   double res;

   // avoid code duplication, therefore only here switch between modified redfield and general foerster
   if(params.modified_redfield==true)
   {
      long_double_complex fac=gppkkpkpk-
                              (gpkpkkpkp-gpkpkkk+long_double(2.0)*long_double_complex(0.0,1.0)*lkpkkpkp)
                              *(gpkpkpkkp-gpkkkkp+long_double(2.0)*long_double_complex(0.0,1.0)*lkpkpkkp);

      long_double_complex exponent=-long_double_complex(0.0,1.0)*(long_double(params.Ek-params.Ekp)/const_hbar+lkkkk+lkpkpkpkp)*long_double(t)
                                   -gkkkk-gkpkpkpkp
                                   +long_double(2.0)*gkpkpkk
                                   +long_double(2.0)*long_double_complex(0.0,1.0)*lkpkpkk*long_double(t);
      res=const_hbar*const_hbar*real(fac*exp(exponent));
   }
   else // generalized foerster
   {
      res=real(exp(-long_double_complex(0.0,1.0)*(long_double(params.Ek-params.Ekp)/const_hbar+lkkkk+lkpkpkpkp)*long_double(t)
                   -gkkkk-gkpkpkpkp));
   }

   // return these two values useful for analytic results for ek==ekp
   params.lkkkk=lkkkk;
   params.re_gpkkkkp=real(gpkkkkp);
   /*
    char s[500];
    sprintf(s,"f%02d%02d.dat",params.ek,params.ekp);
    FILE *fd=fopen(s,"a");
    fprintf(fd,"%e %e %e %e\n",t,res,real(fac),real(exponent));
    fclose(fd);
    */
   return res;
}

double f_combined(double t, void * p)
{
   f_params &params= *reinterpret_cast<f_params *>(p);

   double *lambda;
   double_complex* gamma;

   long_double lkkkk=0.0;
   long_double lkpkpkpkp=0.0;
   long_double lkpkpkk=0.0;
   long_double lkpkkpkp=0.0;
   long_double lkpkpkkp=0.0;

   long_double lkkpkpkp=0.0; // combined

   long_double_complex gkkkk=0.0;
   long_double_complex gkpkpkpkp=0.0;
   long_double_complex gkpkpkk=0.0;
   long_double_complex gpkpkkpkp=0.0;
   long_double_complex gpkpkkk=0.0;
   long_double_complex gpkpkpkkp=0.0;
   long_double_complex gpkkkkp=0.0;

   long_double_complex gpkkpkk=0.0; // combined
   long_double_complex gpkkpkpkp=0.0; // combined

   long_double_complex gppkkpkpk=0.0;

   // sum over all coupled sites
   for(int i=0; i<params.sd_CoupledSite.size(); i++)
   {
      // obtain lambda and gamma for this site
      if(params.sd_Npeaks[i]==1 and params.sd_gammapos_invcm[i][0]==0.0)
      {
         lambda=new double[params.sd_Npeaks[i]];
         gamma=new double_complex[params.sd_Npeaks[i]];
         for(int k=0; k<params.sd_Npeaks[i]; k++)
         {
            lambda[k]=params.sd_lambda_invcm[i][k]*const_invcmtomeV;
            gamma[k]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15), 0.);
         }
      }
      else
      {
         lambda=new double[2*params.sd_Npeaks[i]];
         gamma=new double_complex[2*params.sd_Npeaks[i]];
         for(int k=0; k<params.sd_Npeaks[i]; k++)
         {
            lambda[2*k]=params.sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            lambda[2*k+1]=params.sd_lambda_invcm[i][k]/2.*const_invcmtomeV;
            gamma[2*k]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15),+params.sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
            gamma[2*k+1]=double_complex(+1.0/(params.sd_gamma_fs[i][k]*1.0e-15),-params.sd_gammapos_invcm[i][k]*const_invcmtomeV/const_hbar);
         }
      }

      // compute g,g',g'' at this site at time t (needs to know lambda, gamma, Omega)

      double_complex g,gp,gpp;
      long_double_complex ldg,ldgp,ldgpp;

      ggpgpp_for_sd(t,params.Temp,params.sd_Npeaks[i],lambda,gamma,g,gp,gpp);
      // fprintf(stderr,"t=%e g=(%e,) gp=(%e,) gpp=(%e,)\n",t,real(g),real(gp),real(gpp));
      ldg=g;
      ldgp=gp;
      ldgpp=gpp;

      delete[] lambda;
      delete[] gamma;

      // compute g,g',g'' at this site at time t (needs to know lambda, gamma, Omega)

      int siteindex=params.sd_CoupledSite[i];
      gkkkk    +=((long_double)powl(params.Evk [siteindex],4)*ldg);
      gkpkpkpkp+=((long_double)powl(params.Evkp[siteindex],4)*ldg);
      gkpkpkk  +=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*ldg;

      gpkpkkpkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;
      gpkpkkk  +=long_double(params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*ldgp;
      gpkpkpkkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;
      gpkkkkp  +=long_double(params.Evk [siteindex]*params.Evk [siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*ldgp;

      gpkkpkk  +=long_double(params.Evk [siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*ldgp;
      gpkkpkpkp+=long_double(params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex])*ldgp;

      gppkkpkpk+=long_double(params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex])*ldgpp;

      lkkkk    +=((long_double)powl(params.Evk [siteindex],4)*params.l_list[siteindex]);
      lkpkpkpkp+=((long_double)powl(params.Evkp[siteindex],4)*params.l_list[siteindex]);

      lkpkpkk +=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evk [siteindex])*params.l_list[siteindex];
      lkpkkpkp+=long_double(params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex])*params.l_list[siteindex];
      lkpkpkkp+=long_double(params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evk [siteindex]*params.Evkp[siteindex])*params.l_list[siteindex];

      lkkpkpkp+=long_double(params.Evk [siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex]*params.Evkp[siteindex])*params.l_list[siteindex];
   }
   gkkkk/=const_hbar;

   gkpkpkpkp/=const_hbar;
   gkpkpkk/=const_hbar;
   gpkpkkpkp/=const_hbar;
   gpkpkkk/=const_hbar;
   gpkpkpkkp/=const_hbar;
   gpkkkkp/=const_hbar;

   gpkkpkk  /=const_hbar;
   gpkkpkpkp/=const_hbar;

   gppkkpkpk/=const_hbar;

   lkkkk/=const_hbar;
   lkpkpkpkp/=const_hbar;
   lkpkpkk/=const_hbar;
   lkpkkpkp/=const_hbar;
   lkpkpkkp/=const_hbar;

   lkkpkpkp/=const_hbar;

   double res=0.0;

   {
      long_double_complex fac;
      long_double_complex exponent;

      // part 1
      fac=powl(params.HCW/const_hbar,2)+
          gppkkpkpk-
          (gpkpkkpkp-gpkpkkk+long_double(2.0)*long_double_complex(0.0,1.0)*lkpkkpkp)
          *(gpkpkpkkp-gpkkkkp+long_double(2.0)*long_double_complex(0.0,1.0)*lkpkpkkp);

      exponent=-long_double_complex(0.0,1.0)*(long_double(params.Ek-params.Ekp)/const_hbar+lkkkk+lkpkpkpkp)*long_double(t)
               -gkkkk-gkpkpkpkp
               +long_double(2.0)*gkpkpkk
               +long_double(2.0)*long_double_complex(0.0,1.0)*lkpkpkk*long_double(t);

      res+=const_hbar*const_hbar*real(fac*exp(exponent));

      // part 2
      fac=(params.HCW/const_hbar)*
          ((gpkpkkk-gpkpkkpkp-long_double(2.0)*long_double_complex(0.0,1.0)*lkpkkpkp)
           +(gpkkpkk-gpkkpkpkp-long_double(2.0)*long_double_complex(0.0,1.0)*lkkpkpkp));

      exponent=-long_double_complex(0.0,1.0)*(long_double(params.Ek-params.Ekp)/const_hbar+lkkkk+lkpkpkpkp)*long_double(t)
               -gkkkk-gkpkpkpkp
               +long_double(2.0)*gkpkpkk
               +long_double(2.0)*long_double_complex(0.0,1.0)*lkpkpkk*long_double(t);

      res+=const_hbar*const_hbar*imag(fac*exp(exponent));
   }

   // return these two values useful for analytic results for ek==ekp
   params.lkkkk=lkkkk;
   params.re_gpkkkkp=real(gpkkkkp);
   /*
    char s[500];
    sprintf(s,"f%02d%02d.dat",params.ek,params.ekp);
    FILE *fd=fopen(s,"a");
    fprintf(fd,"%e %e %e %e\n",t,res,real(fac),real(exponent));
    fclose(fd);
    */
   return res;
}

#endif // _HYP2F1_H_

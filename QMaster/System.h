#ifndef SYSTEM_H
#define SYSTEM_H

#include <sstream>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <complex>
#include <vector>
#include "headers.h"

using namespace std;

extern "C" {
   extern void zheev_(char *, char *, int *, double_complex *, int *, double *, double_complex *,
                      int *, double *, int *);
}

// Defines all Systemspecific Parameter.

class System
{
public:
   int Nsites;   //Total number of sites in Hamiltonian
   int Nsites1;  //this is the number of sites EXcluding the |0><0| no exciton ground-state (for pop and coh Nsites=Nsites1, for 1d and 2d Nsites=Nsites1+1
   int Ncoupled; //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
   // define Hamiltonian in site-basis H=(E1+lambda1)|1><1|+(E2+lambda2)|2><2|+J12(|2><1|+|1><2|)
   // see J. Chem. Phys. 130, 234111 (2009), Ishizaki and Fleming
   // define all Parameter
   double invcmtomeV; 	// convert cm^-1 to meV
   double kb;       	// bolzmann constant in meV/K
   double hbar;  		// hbar in meVs
   double Temp; 		// Temperature in K
   double *lambda;         // unit of cm^-1
   Matrix *Ham; 		//Hamiltonian des Exitonsystems in cm^-1

   // calculated in initialize()
   double beta;		// unit 1/meV
   double *Eigenvals;
   double *Eigenvals_full; //eigenvalues full hamiltonian
   Matrix *U;
   Matrix *UT;
   Matrix *U_full;
   Matrix *UT_full;
   Matrix *HamDiag;        //Hamiltonian in Eigenbasis, diagonal HamDiag=U Ham UT
   bool Nan;

   // extra definitions used for 2d spectra
   double_complex *DipoleMoments;
   double_complex *DipoleOperator; // dipole operator
   Matrix *mum, *mup;	// dipole operator exitation and anhilation matrices
   Matrix *mu; 			// dipoole operator mu=mum+mup
private:
   Matrix *Vj;
   Matrix *M1help;
   Matrix *M2help;
   Matrix *M3help;
   Matrix *M4help;
   Matrix *help;
public:
   System(void);
   void initialize(
      int INNsites,   //Total number of Sites in Hamiltonian
      int Ncoupled, //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
      double INinvcmtomeV, 	// convert cm^-1 to meV
      double INkb,       	// bolzmann constant in meV/K
      double INhbar,  		// hbar in meVs
      double INTemp, 		// Temperature in K
      double* INlambdaSum,         // unit of cm^-1
      Matrix* INHam 		//Hamiltonian des Exitonsystems in cm^-1
   );
   void initialize_spectra(
      int INNsites,   //number of excitons (EXCLUDING the |0><0| no-excitation state)
      int Ncoupled, //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
      double INinvcmtomeV, 	// convert cm^-1 to meV
      double INkb,       	// bolzmann constant in meV/K
      double INhbar,  		// hbar in meVs

      double INTemp, 		// Temperature in K
      double* INlambdaSum,         // unit of cm^-1, this is the total sum over all peaks at each site
      Matrix* INHam, 		//Hamiltonian des Exitonsystems in cm^-1

      double_complex* dipoleMoments, // this contains the the dipole moments of all excitons
      double* Pol // laser-light polarization 3-vector
   );
   void computeEigenvects_Eigenvals(void);
   void computeEigenvects_Eigenvals_doubleExciton(void);
   void set_new_laser_direction(double x, double y, double z);
   void create_mum_mup(void);
   void create_mu(void);
   void double_Excitation(void);
   void convertEigenbasis(Matrix* INrho);
   void convertSitebasis(Matrix* INrho);
};

// initialization for calculating coherences and population dynamics
void System::initialize(
   int INNsites,   //Total number of Sites in Hamiltonian
   int INNcoupled, //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
   double INinvcmtomeV, 	// convert cm^-1 to meV
   double INkb,       	// bolzmann constant in meV/K
   double INhbar,  		// hbar in meVs

   double INTemp, 		// Temperature in K
   double* INlambdaSum,         // unit of cm^-1
   Matrix* INHam 		//Hamiltonian des Exitonsystems in cm^-1
)
{
   Nsites=INNsites;
   Nsites1=Nsites; // identical, since no additional |0><0|  involved
   Ncoupled=INNcoupled; //number of sites coupled to phonons times number of Drude-Lorenz peaks per site
   invcmtomeV=INinvcmtomeV;
   kb=INkb;
   hbar=INhbar;
   Temp=INTemp;

   lambda=new double[Nsites];
   for(int i=0; i<Nsites; i++) lambda[i]=INlambdaSum[i];
   Ham=new Matrix(Nsites);
   Ham->assignEntries(*INHam);
   // calculate derived quantities
   beta=1.0/(kb*Temp);

   Eigenvals=new double[Nsites];
   Eigenvals_full=new double[Nsites];
   U=new Matrix(Nsites);
   UT=new Matrix(Nsites);
   U_full=new Matrix(Nsites);
   UT_full=new Matrix(Nsites);
   HamDiag=new Matrix(Nsites);
   help=new Matrix(Nsites);
}



// initialization for calculating 2d spectra
void System::initialize_spectra(
   int INNsites,   //number of excitons (EXCLUDING the |0><0| no-excitation state)
   int INNcoupled, //number of sites coupled to phonon which are evaluated in IF times number of Drude-Lorentzpeaks per site
   double INinvcmtomeV, 	// convert cm^-1 to meV
   double INkb,       	// bolzmann constant in meV/K
   double INhbar,  		// hbar in meVs

   double INTemp, 		// Temperature in K
   double* INlambdaSum,         // unit of cm^-1
   Matrix* INHam, 		//Hamiltonian des Exitonsystems in cm^-1

   double_complex* dipoleMoments, // this contains the the dipole moments of all excitons
   double* Pol // laser-light polarization 3-vector
)
{
   Nsites1=INNsites; // number of excitons (excluding the |0><0| ground state
   Ncoupled=INNcoupled; //number of sites coupled to phonons times number of Drude-Lorenz peaks per site
   invcmtomeV=INinvcmtomeV;
   kb=INkb;
   hbar=INhbar;
   Temp=INTemp;

   beta=1.0/(kb*Temp);

   Nsites=1+Nsites1;  // this is the dimension for the whole system: 0 + 1 Exciton
   // note: the site |0><0| is not counted, since it does not couple

   lambda=new double[Nsites];
   lambda[0]=0.0; // this is for the uncoupled ground state
   for(int i=0; i<Nsites1; i++)
   {
      lambda[i+1]=INlambdaSum[i];
      // attention, shifted by 1 assignement due to prepending |0> state
   }
   DipoleMoments=new double_complex[Nsites1*3]; // store the dipole moments of the sites
   for(int i=0; i<Nsites1; i++)
   {
      DipoleMoments[i*3+0]=dipoleMoments[i*3+0];
      DipoleMoments[i*3+1]=dipoleMoments[i*3+1];
      DipoleMoments[i*3+2]=dipoleMoments[i*3+2];
   }
   DipoleOperator=new double_complex[Nsites];
   DipoleOperator[0]=0.0; // this is for the uncoupled ground state
   for(int i=0; i<Nsites1; i++)
   {
      DipoleOperator[i+1]=Pol[0]*DipoleMoments[i*3+0]+Pol[1]*DipoleMoments[i*3+1]+Pol[2]*DipoleMoments[i*3+2];
      // attention, shifted by 1 assignement due to prepending |0> state
   }

   mum=new Matrix(Nsites);
   mup=new Matrix(Nsites);
   mu= new Matrix (Nsites);
   Ham=new Matrix(Nsites);

   for(int i=0; i<Nsites; i++)
   {
      for(int j=0; j<Nsites; j++)
      {
         if((i==0)||(j==0)) // catches the |0><0| row and column
         {
            Ham->assignEntry(i,j,0.0);
         }
         else
         {
            Ham->assignEntry(i,j,INHam->returnEntry(i-1,j-1));
            // attention: shifted to prepend the |0><0| entries
         }
      }
   }
   // calculate derived quantities
   help=new Matrix(Nsites1); // without |0><0| state
   Eigenvals=new double[Nsites1]; // without |0><0| state
   U=new Matrix(Nsites1); // without |0><0| state
   UT=new Matrix(Nsites1); // without |0><0| state
   HamDiag=new Matrix(Nsites);
   Eigenvals_full=new double[Nsites];
   U_full=new Matrix(Nsites);
   UT_full=new Matrix(Nsites);
}

System::System(void)
{
   Nan=false;
}

void System::computeEigenvects_Eigenvals(void)
{
   char JOBZ='V';
   char UPLO='U';
   int LDA = Nsites1;
   // ANSI C++ requires to init dynamically
   double_complex *WORK;
   WORK=new double_complex[2*Nsites1];
   int LWORK = 2*Nsites1;
   double *RWORK;
   RWORK=new double[3*Nsites1-2];
   int INFO;

   // initialize Eigenvects to Ham which is going to be diagonalized
   // after calling the lapack routine zheev_(...) eigenvects contains all eigenvectors
   int SHIFT=Nsites-Nsites1;
   for(int i=0; i<Nsites1; i++)
   {
      for(int j=0; j<Nsites1; j++)
      {
         U->assignEntry(i,j,Ham->returnEntry(i+SHIFT,j+SHIFT));
      }
   }

   double_complex* doubleU;
   doubleU=new double_complex[Nsites1*Nsites1];

   for(int k=0; k<Nsites1*Nsites1; k++)
   {
      doubleU[k]=(double_complex)(U->entry[k]);
   }
   double *doubleEigenvals;
   doubleEigenvals=new double[Nsites1];
   Eigenvals_full=new double[Nsites];
   for(int k=0; k<Nsites; k++)
   {
      Eigenvals_full[k]=0.;
   }
   for(int k=0; k<Nsites1; k++)
   {
      doubleEigenvals[k]=(double)(Eigenvals[k]);
   }
   // use Lapack to compute eigenvalues and Eigenfunctions, Eigenfunctions are stored in
   // Eigenvects, which was initialized to Ham
   // zheev_ destroyes this initial values
   zheev_(&JOBZ, &UPLO, &Nsites1, doubleU, &LDA, doubleEigenvals, WORK, &LWORK, RWORK, &INFO);
   for(int k=0; k<Nsites1*Nsites1; k++)
   {
      U->entry[k]=(double_complex)(doubleU[k]);
   }
   for(int i=0; i<Nsites1; i++)
   {
      for(int j=0; j<Nsites1; j++)
      {
         U_full->assignEntry(i+SHIFT,j+SHIFT,U->returnEntry(i,j));
      }
   }
   for(int i=0; i<SHIFT; i++)
   {
      U_full->assignEntry(i,i,1.);
   }

   for(int k=0; k<Nsites1; k++)
   {
      Eigenvals[k]=(double)(doubleEigenvals[k]);
      HamDiag->assignEntry(k+SHIFT,k+SHIFT,doubleEigenvals[k]);
      Eigenvals_full[k+SHIFT]=(double)(doubleEigenvals[k]);
   }
   Transpose(UT,U);
   Transpose(UT_full,U_full);
//  cout<<"Exiton Eigenenergies in cm-1"<<endl;
//  for(int i=0; i< Nsites1; i++)
//  {
   //   cout<<Eigenvals[i]/invcmtomeV<<endl;
//  }
//       cout<<"Exiton Eigenenergies in meV"<<endl;
//   for(int i=0; i< Nsites1; i++)
//   {
//        cout<<Eigenvals[i]<<endl;
//   }
   // //   cout<<"Diagonal entries U"<<endl;
//   for(int i=0; i<Nsites1; i++)
//   {
   //     cout<<i<<","<<i<<" :"<<U->returnEntry(i,i)<<endl;
//   }
   //U->print();
   //  MatrixMuliplication(help, Ham, UT, Nsites);
   //  MatrixMuliplication(HamDiag,U,help, Nsites);
   //   cout<<"U Ham UT"<<endl;
   //   HamDiag->print();
   delete[] WORK;
   delete[] RWORK;
}

void System::set_new_laser_direction(double x, double y, double z)
{
   DipoleOperator[0]=0.0; // this is for the uncoupled ground state
   for(int i=0; i<Nsites1; i++)
   {
      DipoleOperator[i+1]=x*DipoleMoments[i*3+0]+y*DipoleMoments[i*3+1]+z*DipoleMoments[i*3+2];
      // attention, shifted by 1 assignement due to prepending |0> state
   }
}

void System::create_mum_mup(void)
{
   for(int i=1; i<=Nsites1; i++)
   {
      cout<<"# dipole moment [site "<<i-1<<"]: "<<DipoleOperator[i]<<endl;
   }

   /// anti-symmetric case:
   for (int i=1; i<=Nsites1; i++)		   // one excited state
   {
      // first the excitations from zero
      mum->assignEntry(0,i,DipoleOperator[i]);
      mup->assignEntry(i,0,DipoleOperator[i]);
   }
}

void System::create_mu(void)
{
//   for(int i=1; i<=Nsites1; i++)
//   {
//      cout<<"DipoleOperator["<<i<<"]: "<<DipoleOperator[i]<<endl;
//   }

   /// anti-symmetric case:
   for (int i=1; i<=Nsites1; i++)		   // one excited state
   {
      // first the excitations from zero
      mu->assignEntry(0,i,DipoleOperator[i]);
      mu->assignEntry(i,0,conj(DipoleOperator[i]));
   }
}

void System::double_Excitation(void)
{
   int li,lj,lk;

   cout<<"# create double excitation"<<endl;
   Nsites=1+Nsites1+Nsites1*(Nsites1-1)/2;  // this is the dimension for the whole system: 0 + 1 Exciton + 2 Excitons

// define the larger lambda
   double *lambdaHelp;
   lambdaHelp=new double[Nsites];
   for(int i=0; i<Nsites; i++)   lambdaHelp[i]=0.0;		// reorganisation energy in meV (bath coupling)
   for(int i=0; i<=Nsites1; i++)   lambdaHelp[i]=lambda[i];		// reorganisation energy in meV (bath coupling)
// assume!!!!!!!!!!!!!!!!!!!!! that the 2 exciton lambdas are simply the sum coupling of double excitation to another bath
   lj=Nsites1+1;  // now lj is in the position of the first 2 exciton states
   for (int i=1; i<Nsites1; i++)
   {
      for (int j=i+1; i<=Nsites1; i++)
      {
         lambdaHelp[lj]=lambda[i]+lambda[j];
         lj++;
      }
   }

// do I need to destroy lambda?
   lambda=new double[Nsites];
   for(int i=0; i<Nsites; i++)   lambda[i]=lambdaHelp[i];

// define the new mum and the new mup:
   mum=new Matrix(Nsites);
   mup=new Matrix(Nsites);
   for (int i=1; i<=Nsites1; i++)		   // one excited state
   {
      // first the excitations from zero
      mum->assignEntry(0,i,DipoleOperator[i]);
      mup->assignEntry(i,0,DipoleOperator[i]);
      // then the excitations from another excited state:
      for (int j=i+1; j<=Nsites1; j++)	 // and the second excited state
      {
         // these matrix elements are at positions (Nsites1+1+ [sum_l=1^(i-1) (Nsites1+1-l)] +j-i)
         lj=Nsites1+j-i;
         for (int l=1; l<i; l++)        lj=lj+Nsites1-l;
         li=Nsites1+i-j;
         for (int l=1; l<j; l++)  li=li+Nsites1-l;

         mum->assignEntry(j,lj,DipoleOperator[i]);
         mum->assignEntry(i,lj,DipoleOperator[j]);    //Boson
         mup->assignEntry(lj,i,DipoleOperator[j]);
         mup->assignEntry(lj,j,DipoleOperator[i]);

//       if (i<j)
//       {
//       	mum->assignEntry(j,lj,-DipoleOperator[i]);
//       	mum->assignEntry(i,lj,DipoleOperator[j]);    //Fermion1
//       	mup->assignEntry(lj,i,DipoleOperator[j]);
//       	mup->assignEntry(lj,j,-DipoleOperator[i]);
//       }
//       if(j<i)
//       {
//         mum->assignEntry(j,lj,DipoleOperator[i]);
//       	mum->assignEntry(i,lj,-DipoleOperator[j]);
//       	mup->assignEntry(lj,i,-DipoleOperator[j]);   //Fermion2
//       	mup->assignEntry(lj,j,DipoleOperator[i]);
//       }
      }
   }


/// and finally H

   Matrix *HamHelp;
   HamHelp=new Matrix(Nsites);
   for (int i=0; i<=Nsites1; i++)
   {
      for (int j=0; j<=Nsites1; j++)	HamHelp->assignEntry(i,j,Ham->returnEntry(i,j));
   }
/// do I need to destroy H?
   Ham=new Matrix(Nsites);
   for (int i=0; i<=Nsites1; i++)
   {
      for (int j=0; j<=Nsites1; j++)	Ham->assignEntry(i,j,HamHelp->returnEntry(i,j));
   }

/// And now add the block for the two exiton Hamiltonian.

//  C++ TRANSLATION from a MATLAB code:
   for (int i=1; i<Nsites1 ; i++)  // one excited state
   {
      for (int j=i+1; j<=Nsites1; j++) // and the second excited state
      {
         // these matrix elements are at positions (Nsites1+1+ [sum_l=1^(i-1) (Nsites1+1-l)] +j-i)
         lj=Nsites1+j-i;
         for (int l=1; l<i; l++)   lj=lj+Nsites1-l;
         // now the diagonal entries are: the matrix elements are: Diagonal ii+jj
         Ham->assignEntry(lj,lj,Ham->returnEntry(j,j)+Ham->returnEntry(i,i));
         // And the off diagonal entries of type i+j and r+k
         for (int r=1; r<=i; r++) // I take an auxillary first exited state
         {
            for (int k=r+1; k<=Nsites1; k++) // and an auxillary second exited state
            {
               lk=Nsites1+k-r;
               for(int l=1; l<r ; l++) 	 lk=lk+Nsites1-l;
               if (i==r && j!=k)
               {
                  Ham->assignEntry(lj,lk,Ham->returnEntry(j,k));
                  Ham->assignEntry(lk,lj,Ham->returnEntry(k,j));
               }
               if (j==k && i!=r)
               {
                  Ham->assignEntry(lk,lj,Ham->returnEntry(r,i));
                  Ham->assignEntry(lj,lk,Ham->returnEntry(i,r));
               }
               if (i==k && j!=r)
               {
                  Ham->assignEntry(lk,lj,Ham->returnEntry(r,j));
                  Ham->assignEntry(lj,lk,Ham->returnEntry(j,r));
               }
            }
         }
      }
   }
//   Ham->print();
}

void System::convertEigenbasis(Matrix* INrho)
{
   Matrix *help;
   help=new Matrix(Nsites);
   MatrixMuliplication(help, INrho, UT_full, Nsites);
   MatrixMuliplication(INrho,U_full, help, Nsites);
}

void System::convertSitebasis(Matrix* INrho)
{
   Matrix *help;
   help=new Matrix(Nsites);
   MatrixMuliplication(help, INrho, U_full, Nsites);
   MatrixMuliplication(INrho,UT_full, help, Nsites);
}



void System::computeEigenvects_Eigenvals_doubleExciton(void)
{
   // rotate along blocks:
   // rotation Matrix of single exciton Manifold+zero exciton state is stored in U_full and UT_full, and eigenvalues in HamDiag
   int Nsites2ex=Nsites1*(Nsites1-1)/2; //sites of double excitons
   Matrix U2ex(Nsites2ex);
   char JOBZ='V';
   char UPLO='U';
   int LDA = Nsites2ex;
   double_complex WORK[2*Nsites2ex];
   int LWORK = 2*Nsites2ex;
   double RWORK[3*Nsites2ex-2];
   int INFO;

   // initialize Eigenvects to Ham which is going to be diagonalized
   // after calling the lapack routine zheev_(...) eigenvects contains all eigenvectors
   int SHIFT=Nsites1+1;
   for(int i=0; i<Nsites2ex; i++)
   {
      for(int j=0; j<Nsites2ex; j++)
      {
         U2ex.assignEntry(i,j,Ham->returnEntry(i+SHIFT,j+SHIFT));
      }
   }

   double_complex* doubleU2ex;
   doubleU2ex=new double_complex[Nsites2ex*Nsites2ex];

   for(int k=0; k<Nsites2ex*Nsites2ex; k++)
   {
      doubleU2ex[k]=(double_complex)(U2ex.entry[k]);
   }
   double *doubleEigenvals2ex;
   doubleEigenvals2ex=new double[Nsites2ex];
   Eigenvals_full=new double[Nsites];  //FIXME if Eigenvals_full needed
   for(int k=0; k<Nsites; k++)
   {
      Eigenvals_full[k]=0.;
   }
   for(int k=0; k<Nsites1+1; k++)
   {
      Eigenvals_full[k]=real(HamDiag->returnEntry(k,k));
   }
   HamDiag=new Matrix(Nsites);
   for(int k=0; k<Nsites2ex; k++)
   {
      doubleEigenvals2ex[k]=0.; //(double)(Eigenvals[k]);
   }
   // use Lapack to compute eigenvalues and Eigenfunctions, Eigenfunctions are stored in
   // Eigenvects, which was initialized to Ham
   // zheev_ destroyes this initial values
   zheev_(&JOBZ, &UPLO, &Nsites2ex, doubleU2ex, &LDA, doubleEigenvals2ex, WORK, &LWORK, RWORK, &INFO);
   for(int k=0; k<Nsites2ex*Nsites2ex; k++)
   {
      U2ex.entry[k]=(double_complex)(doubleU2ex[k]);
   }
   for(int k=0; k<Nsites2ex; k++)
   {
      Eigenvals_full[k+SHIFT]=doubleEigenvals2ex[k];
   }
   U_full=new Matrix(Nsites);
   UT_full=new Matrix(Nsites);
   U_full->assignEntry(0,0,1.);
   for(int i=0; i<Nsites1; i++)
   {
      for(int j=0; j<Nsites1; j++)
      {
         U_full->assignEntry(i+1,j+1,U->returnEntry(i,j));
      }
   }
   for(int i=0; i<Nsites2ex; i++)
   {
      for(int j=0; j<Nsites2ex; j++)
      {
         U_full->assignEntry(i+SHIFT,j+SHIFT,U2ex.returnEntry(i,j));
      }
   }
   for(int k=0; k<Nsites; k++)
   {
      HamDiag->assignEntry(k,k,Eigenvals_full[k]);
//      cout<<HamDiag->returnEntry(k,k)<<endl;
   }

   Transpose(UT_full,U_full);
}


#endif

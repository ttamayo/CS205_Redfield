#ifndef LIOUVILLEPHONON_SECULARREDFIELD_H
#define LIOUVILLEPHONON_SECULARREDFIELD_H

//#include<gsl/gsl_complex.h>
//#include<gsl/gsl_complex_math.h>
//#include<gsl/gsl_sf.h>

#include "headers.h"

class LiouvillePhonon_secularRedfield: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
private:
   vector<double> listlambda;
   vector<double_complex> listgamma;
   vector<int> listNDL;
   vector<int> listNDLsite;
   vector<Matrix*> Gamma_Rate;
   Matrix *help1;
   Matrix *rho_new;
   Matrix *rho_old;
   myreal *h_rhoEntries_new;
   myreal *h_rhoEntries_old;
public:
   LiouvillePhonon_secularRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo);
   ~LiouvillePhonon_secularRedfield();
   void AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma); //INNDL stands for the total number of Drude-Lorentz peaks in the spectral density for site INsite
   void initialize();  //init Operator Lambda and Lambdadag
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory();
private:
   double phononStatistic(double energy);
   double LDshifted(double lambda, double gamma, double shiftgamma, double omega);
   double SpectralDens(int j, double energy);
   double dSpectralDens(int j);   // dJ(E)/dE|_E=0
   double gamma(int j, double energy);
   void createGamma(int j, Matrix *INGamma_Rate);
   double sqabs(double_complex c); //computes |c|^2
   void  D(int m, int M, int N, Matrix *INrho, Matrix *result);
};

LiouvillePhonon_secularRedfield::LiouvillePhonon_secularRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo)
{
   help1=new Matrix(system->Nsites);
   rho_new=new Matrix(system->Nsites);
   rho_old=new Matrix(system->Nsites);
   h_rhoEntries_new=new double[2*system->Nsites*system->Nsites];
   h_rhoEntries_old=new double[2*system->Nsites*system->Nsites];
}

LiouvillePhonon_secularRedfield::~LiouvillePhonon_secularRedfield()
{
}

void LiouvillePhonon_secularRedfield::AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma)
{
   listNDL.push_back(INNDL);
   listNDLsite.push_back(INsite);
   for(int i=0; i<INNDL; i++)
   {
      listlambda.push_back(INlambda[i]/2.); //in this routine the spectral density is normalized to 2*lambda. Thus 1/2
      listgamma.push_back(INgamma[i]);
//   cout<<"site "<<INsite<<" lambda "<<INlambda[i]/system->invcmtomeV<<" invgamma in fs "<<1.e15/INgamma[i]<<endl;
   }
   Matrix* Gam;
   Gam=new Matrix(system->Nsites);
   Gamma_Rate.push_back(Gam);
}


double LiouvillePhonon_secularRedfield::phononStatistic(double energy)
{
   double mu=0;
   return 1./(exp(system->beta*(energy-mu))-1);
}

double LiouvillePhonon_secularRedfield::LDshifted(double lambda, double gamma, double shift, double omega)
{
   //in this routine the spectral density is normalized to 2*lambda. To be consistent, the input lambdas in AddSite were multiplied by 1/2, see above
   double J;
   J=2*lambda*(omega*gamma)/((omega-shift)*(omega-shift)+gamma*gamma)/system->hbar;
   J+=2*lambda*(omega*gamma)/((omega+shift)*(omega+shift)+gamma*gamma)/system->hbar;
   return J;
}

double LiouvillePhonon_secularRedfield::SpectralDens(int j, double energy)   //FIXME what is j, not sites should come from listNDLsite
{
   double omega=energy/system->hbar;
   if(energy>0)
   {
      double J=0.;
      int startpos=0;
      for(int k=0; k<j; k++)
      {
         startpos+=listNDL[k];
      }
      for(int i=0; i<listNDL[j]; i++)
      {
         J+=LDshifted(listlambda[startpos+i], real(listgamma[startpos+i]) , imag(listgamma[startpos+i]), omega);
      }
      return J;
   }
   else
   {
      return 0.;
   }
}

double LiouvillePhonon_secularRedfield::dSpectralDens(int j)
{

   double dJ=0.;
   int startpos=0;
   for(int k=0; k<j; k++)
   {
      startpos+=listNDL[k];
   }
   for(int i=0; i<listNDL[j]; i++)
   {
      dJ+=4.*listlambda[startpos+i]*real(listgamma[startpos+i])/(real(listgamma[startpos+i])*real(listgamma[startpos+i])+imag(listgamma[startpos+i])*imag(listgamma[startpos+i]));
   }
   return dJ/system->hbar/system->hbar;
}

double LiouvillePhonon_secularRedfield::gamma(int j, double energy)
{
   double Repart=0;
   if(energy>0)
   {
      Repart=2*M_PI*SpectralDens(j,energy)*(phononStatistic(energy)); //FIXME 2*PI?
   }
   else if(energy<0)
   {
      //FIXME 2*PI? only factor 2* needed, factor PI is already in spectral densiy,
      // note different definition of spectral density of Non-Markovian Quantum Jumps in Excitonic Energy Transfer eq.(11)
      // and definition in Quantum transport through complex networks - from light-harvesting proteins to semiconductor devices eq.(3.5)
      Repart=2*M_PI*SpectralDens(j,-energy)*(phononStatistic(-energy)+1); // Note phononStatistic(energy)+1=-phononStatistic(-energy)
   }
   else if(energy==0)
   {
      //FIXME 2*PI?
      Repart=2*M_PI*system->kb*system->Temp*dSpectralDens(j);
   }
   return Repart;
}

// geht ueber in ein Create rate
void LiouvillePhonon_secularRedfield::createGamma(int j,  Matrix *INGamma_Rate)
{
   double rate;
   for(int M=0; M<system->Nsites; M++)  //Loop over eigenstates |M>
   {
      for(int N=0; N<system->Nsites; N++) //Loop over eingenstates <N|
      {
         rate=gamma(j,(system->Eigenvals_full[M]-system->Eigenvals_full[N])); //FIXME this needs to be adaped,
         INGamma_Rate->assignEntry(M,N,double_complex(rate,0.));
      }
   }
}

void LiouvillePhonon_secularRedfield::initialize()
{
   for(int i=0; i<listNDLsite.size(); i++)
   {
      createGamma(i,Gamma_Rate[i]);
   }
}

double LiouvillePhonon_secularRedfield::sqabs(double_complex c)
{
   return real(c)*real(c)+imag(c)*imag(c);
}

///computes result=D[A_m(EM-EN)]*INrho
void LiouvillePhonon_secularRedfield::D(int m, int M, int N, Matrix *INrho, Matrix *result)
{
   double_complex cmM;
   double_complex cmN;
   cmM=sqabs(system->U_full->returnEntry(M,m));
   cmN=sqabs(system->U_full->returnEntry(N,m));
//
   result->assignEntry(M, M, cmN*cmM*INrho->returnEntry(N,N)); //A_m(EM-EN)rho A_m(EM-EN)^dag

   for(int k=0; k<system->Nsites; k++)                             //A_j(E_M) A_j^dag(E_M) rho
   {
      result->EntryAdd(N, k, -0.5*cmM*cmN*INrho->returnEntry(N,k));
   }
   for(int k=0; k<system->Nsites; k++)                             //rho A_j(E_M) A_j^dag(E_M)
   {
      result->EntryAdd(k, N, -0.5*cmM*cmN*INrho->returnEntry(k,N));
   }
}

/// Liouville Operator, Exciton-Photon Coupling in secular Redfield
void LiouvillePhonon_secularRedfield::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   // load from device memory to CPU memory (secular Redfield part is evaluated on CPU, not-parallel)
   clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);
   clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmaold), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_old,0, NULL, NULL);
   for(int i=0; i<system->Nsites*system->Nsites; i++)
   {
      rho_new->entry[i]=double_complex(h_rhoEntries_new[2*i],h_rhoEntries_new[2*i+1]);
      rho_old->entry[i]=double_complex(h_rhoEntries_old[2*i],h_rhoEntries_old[2*i+1]);
   }
   for(int m=0; m<listNDLsite.size(); m++) // loop over all coupled sites
   {
      for(int M=0; M<system->Nsites; M++)  //Loop over eigenstates |M>
      {
         for(int N=0; N<system->Nsites; N++) //Loop over eingenstates <N|
         {
            D(listNDLsite[m], M, N, rho_old, help1);
            double_complex rate;
            /// FIXME define Matrix-vector length m,  gammaPhonon(m,M,N)
            rate=dt*Gamma_Rate[m]->returnEntry(M,N);
            MultiplyScalar(help1,rate, system->Nsites);
            MatrixAdd(rho_new, help1, system->Nsites);
            help1->clearEntries();
         }
      }
   }
   // memcopy of updated new density matrix from CPU to Device
   for(int i=0; i<system->Nsites*system->Nsites; i++)
   {
      h_rhoEntries_new[2*i]=real(rho_new->entry[i]);
      h_rhoEntries_new[2*i+1]=imag(rho_new->entry[i]);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);
}

//

void LiouvillePhonon_secularRedfield::freeDeviceMemory()
{
//	   listlambda.clear();
//	   listgamma.clear();
//	   listNDL.clear();
//	   listNDLsite.clear();
//	   Lambda.clear();
//	   Lambdadag.clear();
//	   Vjdiag.clear();
//	   delete[] help1;
//	   delete[] help2;
//	   delete[] help3;
//	   delete[] help4;
//	   delete[] help5;
//	   delete[] rho_new;
//	   delete[] rho_old;
//	   delete[] h_rhoEntries_new;
//	   delete[] h_rhoEntries_old;
}


#endif

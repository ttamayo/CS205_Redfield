#ifndef LIOUVILLEPHONONDOUBLEEXCITON_FULLREDFIELD_H
#define LIOUVILLEPHONONDOUBLEEXCITON_FULLREDFIELD_H

#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_sf.h>

#include "headers.h"

class LiouvillePhononDouble_fullRedfield: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
private:
   vector<double> listlambda;
   vector<double_complex> listgamma;
   vector<int> listNDL;
   vector<int> listNDLsite;
   vector<Matrix*> Lambda; //
   vector<Matrix*> Lambdadag;
   vector<Matrix*>	Vjdiag;
   Matrix *help1;
   Matrix *help2;
   Matrix *help3;
   Matrix *help4;
   Matrix *help5;
   Matrix *rho_new;
   Matrix *rho_old;
   myreal *h_rhoEntries_new;
   myreal *h_rhoEntries_old;
   myreal hbardt;
public:
   LiouvillePhononDouble_fullRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo);
   ~LiouvillePhononDouble_fullRedfield();
   void AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma); //INNDL stands for the total number of Drude-Lorentz peaks in the spectral density for site INsite
   void initialize();  //init Operator Lambda and Lambdadag
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory();
private:
   double_complex Digamma(double_complex z);
   double phononStatistic(double energy);
   double LDshifted(double lambda, double gamma, double shiftgamma, double omega);
   double SpectralDens(int j, double energy);
   double dSpectralDens(int j);   // dJ(E)/dE|_E=0
   double_complex gamma(int j, double energy);
   void createLambda(int j, Matrix *INLambda, Matrix *INLambdadag, Matrix* INVjdiag);
   void Liouv(Matrix *INrhoalt, Matrix *INrhonew, myreal dt);
   int id2ex(int k, int j, int Nsites1);
};


int LiouvillePhononDouble_fullRedfield::id2ex(int k, int j, int Nsites1) //returns index in Hamiltonian for double excitation |j,k> if j>k
{
   int id=Nsites1+(j-k);
   for(int i=1; i<k; i++)
   {
      id+=Nsites1-i;
   }
   return id;
}

LiouvillePhononDouble_fullRedfield::LiouvillePhononDouble_fullRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo)
{
   help1=new Matrix(system->Nsites);
   help2=new Matrix(system->Nsites);
   help3=new Matrix(system->Nsites);
   help4=new Matrix(system->Nsites);
   help5=new Matrix(system->Nsites);
   rho_new=new Matrix(system->Nsites);
   rho_old=new Matrix(system->Nsites);
   h_rhoEntries_new=new double[2*system->Nsites*system->Nsites];
   h_rhoEntries_old=new double[2*system->Nsites*system->Nsites];
}

LiouvillePhononDouble_fullRedfield::~LiouvillePhononDouble_fullRedfield()
{
}



void LiouvillePhononDouble_fullRedfield::AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma)
{
   listNDL.push_back(INNDL);
   listNDLsite.push_back(INsite);
   for(int i=0; i<INNDL; i++)
   {
      listlambda.push_back(INlambda[i]/2.); //in this routine the spectral density is normalized to 2*lambda. Thus 1/2
      listgamma.push_back(INgamma[i]);
//   cout<<"site "<<INsite<<" lambda "<<INlambda[i]/system->invcmtomeV<<" invgamma in fs "<<1.e15/INgamma[i]<<endl;
   }
   Matrix* L;
   L=new Matrix(system->Nsites);
   Matrix* Ldag;
   Ldag=new Matrix(system->Nsites);
   Lambda.push_back(L);
   Lambdadag.push_back(Ldag);
   Matrix* Vj;
   Vj=new Matrix(system->Nsites);
   Vjdiag.push_back(Vj);
}

double_complex LiouvillePhononDouble_fullRedfield::Digamma(double_complex z)
{
   double x,y,x0,x1,y1,psr,psi;
   double q0,q2,rr,ri,th,tn,tm,ct2;
   double_complex ps;

   int n,k;
   static double a[] =
   {
      -0.8333333333333e-01,
      0.83333333333333333e-02,
      -0.39682539682539683e-02,
      0.41666666666666667e-02,
      -0.75757575757575758e-02,
      0.21092796092796093e-01,
      -0.83333333333333333e-01,
      0.4432598039215686
   };

   x = real(z);
   y = imag(z);
   x1 = 0.0;
   n = 0;
   if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
   {
      ps = double_complex (1e308,0);
   }
   else
   {
      if (x < 0.0)
      {
         x1 = x;
         y1 = y;
         x = -x;
         y = -y;
      }
      x0 = x;
      if (x < 8.0)
      {
         n = 8 - (int)x;
         x0 = x + n;
      }
      if ((x0 == 0.0) && (y != 0.0))
         th = 0.5*M_PI;
      if (x0 != 0.0)
         th = atan(y/x0);
      q2 = x0*x0+y*y;
      q0 = sqrt(q2);
      psr = log(q0)-0.5*x0/q2;
      psi = th+0.5*y/q2;
      for (k=1; k<=8; k++)
      {
         psr += (a[k-1]*pow(q2,-k)*cos(2.0*k*th));
         psi -= (a[k-1]*pow(q2,-k)*sin(2.0*k*th));
      }
      if (x < 8.0)
      {
         rr = 0.0;
         ri = 0.0;
         for (k=1; k<=n; k++)
         {
            rr += ((x0-k)/(pow(x0-k,2.0)+y*y));
            ri += (y/(pow(x0-k,2.0)+y*y));
         }
         psr -= rr;
         psi += ri;
      }
      if (x1 < 0.0)
      {
         tn = tan(M_PI*x);
         tm = tanh(M_PI*y);
         ct2 = tn*tn+tm*tm;
         psr = psr+x/(x*x+y*y)+M_PI*(tn-tn*tm*tm)/ct2;
         psi = psi-y/(x*x+y*y)-M_PI*tm*(1.0+tn*tn)/ct2;
         x = x1;
         y = y1;
      }
      ps = double_complex(psr,psi);
   }
   return ps;
}

double LiouvillePhononDouble_fullRedfield::phononStatistic(double energy)
{
   double mu=0;
   return 1./(exp(system->beta*(energy-mu))-1);
}

double LiouvillePhononDouble_fullRedfield::LDshifted(double lambda, double gamma, double shift, double omega)
{
   //in this routine the spectral density is normalized to 2*lambda. To be consistent, the input lambdas in AddSite were multiplied by 1/2, see above
   double J;
   J=2*lambda*(omega*gamma)/((omega-shift)*(omega-shift)+gamma*gamma)/system->hbar;
   J+=2*lambda*(omega*gamma)/((omega+shift)*(omega+shift)+gamma*gamma)/system->hbar;
   return J;
}

double LiouvillePhononDouble_fullRedfield::SpectralDens(int j, double energy)   //FIXME what is j, not sites should come from listNDLsite
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

double LiouvillePhononDouble_fullRedfield::dSpectralDens(int j)
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

double_complex LiouvillePhononDouble_fullRedfield::gamma(int j, double energy)
{
   double hbar=system->hbar;
   double omega=energy/hbar;
   double beta=system->beta;
   if(energy>0)
   {
      double Repart=0.;
      double Impart=0.;
      int startpos=0;
      for(int k=0; k<j; k++)
      {
         startpos+=listNDL[k];
      }
      for(int i=0; i<listNDL[j]; i++)
      {
         double ImpartDLpeak=0;
         double gamma=real(listgamma[startpos+i]);
         double lambda=listlambda[startpos+i];
         double Omega=imag(listgamma[startpos+i]);

         double square=gamma*gamma+omega*omega+Omega*Omega;

         ImpartDLpeak=M_PI*gamma/(2.*omega);
         ImpartDLpeak+=M_PI/(beta*hbar*gamma);
         ImpartDLpeak+=-2.*M_PI*Omega*Omega/(hbar*beta*gamma*square);
         ImpartDLpeak+=M_PI*gamma*Omega*Omega*(1.-(omega*omega-Omega*Omega)/(gamma*gamma))/(2.*omega*square);
         double_complex argDigamma1=double_complex(hbar*beta*gamma/(2.*M_PI),-hbar*beta*Omega/(2.*M_PI));
         double_complex psi1=Digamma(argDigamma1);
         ImpartDLpeak+=-Omega/gamma*(gamma*gamma+Omega*Omega-omega*omega)/square*imag(psi1)+real(psi1);
         ImpartDLpeak+=-gsl_sf_psi_1piy(beta*hbar*omega/(2.*M_PI));
         ImpartDLpeak*=-LDshifted(lambda, gamma, Omega, energy/system->hbar)/(M_PI);
         Impart+=ImpartDLpeak;
      }
      Repart=SpectralDens(j,energy)*(phononStatistic(energy)+1);
      return double_complex(Repart,Impart);
   }
   else if(energy<0)
   {
      double Repart=0.;
      double Impart=0.;
      int startpos=0;
      for(int k=0; k<j; k++)
      {
         startpos+=listNDL[k];
      }
      for(int i=0; i<listNDL[j]; i++)
      {
         double ImpartDLpeak=0;
         double gamma=real(listgamma[startpos+i]);
         double lambda=listlambda[startpos+i];
         double Omega=imag(listgamma[startpos+i]);

         double square=gamma*gamma+omega*omega+Omega*Omega;

         ImpartDLpeak=M_PI*gamma/(2.*omega);
         ImpartDLpeak+=M_PI/(beta*hbar*gamma);
         ImpartDLpeak+=-2.*M_PI*Omega*Omega/(hbar*beta*gamma*square);
         ImpartDLpeak+=M_PI*gamma*Omega*Omega*(1.-(omega*omega-Omega*Omega)/(gamma*gamma))/(2.*omega*square);
         double_complex argDigamma1=double_complex(hbar*beta*gamma/(2.*M_PI),-hbar*beta*Omega/(2.*M_PI));
         double_complex psi1=Digamma(argDigamma1);
         ImpartDLpeak+=-Omega/gamma*(gamma*gamma+Omega*Omega-omega*omega)/square*imag(psi1);
         ImpartDLpeak+=real(psi1);
         ImpartDLpeak+=-gsl_sf_psi_1piy(beta*hbar*omega/(2.*M_PI));
         ImpartDLpeak*=LDshifted(lambda, gamma, Omega, -energy/system->hbar)/(M_PI);
         Impart+=ImpartDLpeak;
      }
//   Impart*=SpectralDens(j, -energy)/(M_PI);
      Repart=SpectralDens(j,-energy)*(phononStatistic(-energy)); //Note phononStatistic(energy)+1=-phononStatistic(-energy)
      return double_complex(Repart,Impart);
   }
   else if(energy==0)
   {
      double Repart=system->kb*system->Temp*dSpectralDens(j);
      double Impart=0.;
      int startpos=0;
      for(int k=0; k<j; k++)
      {
         startpos+=listNDL[k];
      }
      for(int i=0; i<listNDL[j]; i++)
      {
         double lambda=listlambda[startpos+i];
         Impart+=-2.*lambda/hbar;
      }
      return double_complex(Repart, Impart);
   }
}

void LiouvillePhononDouble_fullRedfield::createLambda(int j, Matrix *INLambda, Matrix *INLambdadag, Matrix* INVjdiag)
{
   int i=0;
   for(int M=0; M<system->Nsites; M++)
   {
      for(int N=0; N<system->Nsites; N++)
      {

         double_complex c;
         double_complex cconj;
         double_complex cVjdiag;
         double_complex entryVj=0;
         double_complex entryLambda=0;
         double_complex entryLambdadag=0;
         c=gamma(j,system->Eigenvals_full[N]-system->Eigenvals_full[M])*system->UT_full->returnEntry(listNDLsite[j],M)*system->U_full->returnEntry(N,listNDLsite[j]);
         cconj=conj(c);
         cVjdiag=conj(system->U_full->returnEntry(M,listNDLsite[j]))*system->U_full->returnEntry(N, listNDLsite[j]);

         entryLambda+=c;
         entryLambdadag+=cconj;
         entryVj+=cVjdiag;
         for(int kk=1; kk<=system->Nsites1; kk++)
         {
            if(listNDLsite[j]<kk)
            {
               int idex=id2ex(listNDLsite[j], kk, system->Nsites1);
               c=gamma(j,system->Eigenvals_full[N]-system->Eigenvals_full[M])*system->UT_full->returnEntry(idex,M)*system->U_full->returnEntry(N,idex);
               cconj=conj(c);
               cVjdiag=conj(system->U_full->returnEntry(M,idex))*system->U_full->returnEntry(N, idex);
               entryLambda+=c;
               entryLambdadag+=cconj;
               entryVj+=cVjdiag;
            }
            if(listNDLsite[j]>kk)
            {
               int idex=id2ex(kk,listNDLsite[j], system->Nsites1);
               c=gamma(j,system->Eigenvals_full[N]-system->Eigenvals_full[M])*system->UT_full->returnEntry(idex,M)*system->U_full->returnEntry(N,idex);
               cconj=conj(c);
               cVjdiag=conj(system->UT_full->returnEntry(idex,M))*system->U_full->returnEntry(N, idex);
               entryLambda+=c;
               entryLambdadag+=cconj;
               entryVj+=cVjdiag;
            }
         }
         INLambda->assignEntry(M,N, entryLambda);
         INLambdadag->assignEntry(N,M, entryLambdadag);
         INVjdiag->assignEntry(M,N,entryVj);

      }
   }
}

void LiouvillePhononDouble_fullRedfield::initialize()
{
   for(int i=0; i<listNDLsite.size(); i++)
   {
      createLambda(i,Lambda[i],Lambdadag[i], Vjdiag[i]);
//  cout<<"# Lambdadag["<<i<<"]"<<endl;
//  Lambdadag[i]->print();
//  cout<<"# Vjdiag["<<i<<"]"<<endl;
//  Vjdiag[i]->print();
   }
}

void LiouvillePhononDouble_fullRedfield::Liouv(Matrix *INrhoalt, Matrix *INrhonew, myreal dt)
{
   for(int i=0; i<listNDLsite.size(); i++) // loop over all couplings
   {
      //compute help2=Vj Lambda rhodiag
      MatrixMuliplication(help1, Lambda[i], INrhoalt, system->Nsites);
      MatrixMuliplication(help2, Vjdiag[i] , help1, system->Nsites);
      //compute help3=Vj*rhodiag*Lambdadag
      MatrixMuliplication(help1, INrhoalt, Lambdadag[i], system->Nsites);
      MatrixMuliplication(help3, Vjdiag[i] , help1, system->Nsites);
      //compute help4=Lambda*rhodiag*Vj
      MatrixMuliplication(help1, INrhoalt, Vjdiag[i], system->Nsites);
      MatrixMuliplication(help4, Lambda[i], help1, system->Nsites);
//   compute help5=rhodiag*Lambdadag*Vj
      MatrixMuliplication(help1, Lambdadag[i], Vjdiag[i], system->Nsites);
      MatrixMuliplication(help5, INrhoalt, help1, system->Nsites);
      //compute help3=help3+help4-help5-help2
      MatrixAdd(help3, help4, system->Nsites);
      MatrixSubtract(help3, help5, system->Nsites);
      MatrixSubtract(help3, help2, system->Nsites);
      //Add help3 to INrhonew
      MultiplyScalar(help3,dt,system->Nsites);
      MatrixAdd(INrhonew,help3, system->Nsites);
   }
}

void LiouvillePhononDouble_fullRedfield::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);
   clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmaold), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_old,0, NULL, NULL);
   for(int i=0; i<system->Nsites*system->Nsites; i++)
   {
      rho_new->entry[i]=double_complex(h_rhoEntries_new[2*i],h_rhoEntries_new[2*i+1]);
      rho_old->entry[i]=double_complex(h_rhoEntries_old[2*i],h_rhoEntries_old[2*i+1]);
   }
   Liouv(rho_old,rho_new,dt);
   for(int i=0; i<system->Nsites*system->Nsites; i++)
   {
      h_rhoEntries_new[2*i]=real(rho_new->entry[i]);
      h_rhoEntries_new[2*i+1]=imag(rho_new->entry[i]);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);

}

void LiouvillePhononDouble_fullRedfield::freeDeviceMemory()
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

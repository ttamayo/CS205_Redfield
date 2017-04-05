#ifndef GETDENSITYMATRIXSITEBASISDIAG_H
#define GETDENSITYMATRIXSITEBASISDIAG_H


#include "headers.h"

class GetDensityMatrixSiteBasisDiag: public Observable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   int nsteps;
   vector<double_complex> data;
   int Ndata;
private:
   Matrix *rho;
   double *h_rhoEntries;
public:
   GetDensityMatrixSiteBasisDiag(System &INsystem,OpenCL_Init INOpenCLinfo, int INt_mod, int INNsteps);
   void execute(double time, int it,  cl_mem d_sigma_System);
};

GetDensityMatrixSiteBasisDiag::GetDensityMatrixSiteBasisDiag(System &INsystem, OpenCL_Init INOpenCLinfo,int INt_mod, int INNsteps) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), t_mod(INt_mod), nsteps(INNsteps)
{
   rho=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
   Ndata=0;
}


void GetDensityMatrixSiteBasisDiag::execute(double time, int it,  cl_mem d_sigma_System)
{
   if(it%t_mod==0)
   {
      Ndata+=1;
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=double_complex(h_rhoEntries[2*i],h_rhoEntries[2*i+1]);
      }
      for(int i=0; i<system->Nsites; i++)
      {
         data.push_back(rho->returnEntry(i,i));
         //check if nan, if yes no further output, program aborts
         double test=abs(data.back());
         if(test != test)
         {
            cout<<"overflow: "<<test<<endl;
            system->Nan=true;
         }
      }
   }
}

#endif

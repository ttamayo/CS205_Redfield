#ifndef GetDensityMatrixElement_EigenBasis_Redfield_REDFIELD_H
#define GetDensityMatrixElement_EigenBasis_Redfield_REDFIELD_H


#include "headers.h"

/// Returns Element |Ei >< Ej| of density matrix in Eigenbasis
/// Achtung bei Numerik Energieen gehen von E0 bis E6, bei Benennung der Filenamenen
/// addieren wir +1 so dass Energieen von E1 bis E7
class GetDensityMatrixElement_EigenBasis_Redfield: public Observable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   int Ei;
   int Ej;
   vector<double_complex> data;
private:
   Matrix *rho;
   Matrix *rhoEigBasis;
   Matrix *Help1;
   double *h_rhoEntries;
public:
   GetDensityMatrixElement_EigenBasis_Redfield(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INEi, int INEj);
   ~GetDensityMatrixElement_EigenBasis_Redfield(void);
   void execute(double time, int it, cl_mem d_sigma_System);
};

GetDensityMatrixElement_EigenBasis_Redfield::GetDensityMatrixElement_EigenBasis_Redfield(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INEi, int INEj) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), t_mod(INt_mod), Ei(INEi), Ej(INEj)
{
   rho=new Matrix(system->Nsites);
   rhoEigBasis=new Matrix(system->Nsites);
   Help1=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
}

GetDensityMatrixElement_EigenBasis_Redfield::~GetDensityMatrixElement_EigenBasis_Redfield(void)
{
   free(rho);
   free(rhoEigBasis);
   free(Help1);
   free(h_rhoEntries);
}

void GetDensityMatrixElement_EigenBasis_Redfield::execute(double time, int it, cl_mem d_sigma_System)
{
   if(it%t_mod==0)
   {
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=h_rhoEntries[2*i]+double_complex(0.,h_rhoEntries[2*i+1]);
      }
      //
      data.push_back(rho->returnEntry(Ei, Ej));
      double test=abs(data.back());
      if(test != test)
      {
         cout<<"overflow: "<<test<<endl;
         system->Nan=true;
      }
      //
   }
}

#endif


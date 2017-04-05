#ifndef GETDENSITYMATRIXELEMENT_EIGENBASIS_H
#define GETDENSITYMATRIXELEMENT_EIGENBASIS_H


#include "headers.h"

/// Returns Element |Ei >< Ej| of density matrix in Eigenbasis
/// Achtung bei Numerik Energieen gehen von E0 bis E6, bei Benennung der Filenamenen
/// addieren wir +1 so dass Energieen von E1 bis E7
class GetDensityMatrixElement_EigenBasis: public Observable
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
   bool PrintAllDiagEntries;
public:
   GetDensityMatrixElement_EigenBasis(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INEi, int INEj, bool INPrintAllDiagEntries);
   ~GetDensityMatrixElement_EigenBasis(void);
   void execute(double time, int it, cl_mem d_sigma_System);
};

GetDensityMatrixElement_EigenBasis::GetDensityMatrixElement_EigenBasis(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INEi, int INEj, bool INPrintAllDiagEntries) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), t_mod(INt_mod), Ei(INEi), Ej(INEj), PrintAllDiagEntries(INPrintAllDiagEntries)
{
   rho=new Matrix(system->Nsites);
   rhoEigBasis=new Matrix(system->Nsites);
   Help1=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
}

GetDensityMatrixElement_EigenBasis::~GetDensityMatrixElement_EigenBasis(void)
{
   free(rho);
   free(rhoEigBasis);
   free(Help1);
   free(h_rhoEntries);
}

void GetDensityMatrixElement_EigenBasis::execute(double time, int it, cl_mem d_sigma_System)
{
   if(it%t_mod==0)
   {
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=h_rhoEntries[2*i]+double_complex(0.,h_rhoEntries[2*i+1]);
      }
      //
      MatrixMuliplication(Help1, rho, system->UT, system->Nsites);
      MatrixMuliplication(rhoEigBasis,system->U,Help1, system->Nsites);

      if(PrintAllDiagEntries==false)
      {
         data.push_back(rhoEigBasis->returnEntry(Ei, Ej));
         double test=abs(data.back());
         if(test != test)
         {
            cout<<"overflow: "<<test<<endl;
            system->Nan=true;
         }
      }
      if(PrintAllDiagEntries==true)
      {
         for(int i=0; i<system->Nsites; i++)
         {
            for(int j=0; j<system->Nsites; j++)
            {
               data.push_back(rhoEigBasis->returnEntry(i,j));
               double test=abs(data.back());
               if(test != test)
               {
                  cout<<"overflow: "<<test<<endl;
                  system->Nan=true;
               }
            }
         }
      }
   }
}

#endif

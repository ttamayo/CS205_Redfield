#ifndef BREAKBYNORM_ANTIHERM_H
#define BREAKBYNORM_ANTIHERM_H

#include "headers.h"

class BreakByNorm_antiHerm: public BreakRule
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   int breakID; // id of break condition, note that there might be several break conditions
private:
   Matrix *rho;
   double *h_rhoEntries;
public:
   BreakByNorm_antiHerm(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INbreakID);
   void execute(double *breakval, cl_mem d_sigma_System, int it);
};

BreakByNorm_antiHerm::BreakByNorm_antiHerm(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, int INbreakID) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo),  t_mod(INt_mod), breakID(INbreakID)
{
   rho=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
}


void BreakByNorm_antiHerm::execute(double *breakval, cl_mem d_sigma_System, int it)
{
   if(it%t_mod==0)
   {
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=h_rhoEntries[2*i]+double_complex(0.,h_rhoEntries[2*i+1]);
      }
      double norm=0.;
      for(int i=0; i<system->Nsites; i++)
      {
         norm+=real(rho->returnEntry(i, i));
      }
      breakval[breakID]=norm;
// cout<<"set breakval"<<endl;
      cout<<"# norm="<<scientific<<breakval[breakID]<<endl;
      cout.flush();
   }
}

#endif

#ifndef RETURNTRACE_MUMRHO_H
#define RETURNTRACE_MUMRHO_H

#include "headers.h"

// returns real and imaginary part of tr(mum rho(t)) as vectors

class ReturnTrace_mumRho: public ManipulatingObservable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   Matrix *mum;
   Matrix *Help;
   int delayTime;
   FILE *OUTr;
   FILE *OUTi;
   FILE *OUT_dat;
   vector<double> tracer;
   vector<double> tracei;
private:
   Matrix *rho;
   double *h_rhoEntries;
public:

   ReturnTrace_mumRho(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, Matrix *INmum, Matrix *INHelp, int INdelayTime); // returns the trace of the time evolution to the main programme
   void execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help);
};

ReturnTrace_mumRho::ReturnTrace_mumRho(System &INsystem, OpenCL_Init INOpenCLinfo, int INt_mod, Matrix *INmum, Matrix *INHelp, int INdelayTime) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), t_mod(INt_mod), mum(INmum), Help(INHelp), delayTime(INdelayTime)
{
   rho=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
}


void ReturnTrace_mumRho::execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help)
{

   if((it>=delayTime)&&((it-delayTime)%t_mod==0))
   {
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(double)*2*system->Nsites*system->Nsites, h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=double_complex(h_rhoEntries[2*i],h_rhoEntries[2*i+1]);
      }
      MatrixMuliplication(Help, mum, rho, system->Nsites); // first=second * third

      double Repart=0.0;
      double Impart=0.0;
      for (int i=0; i<system->Nsites; i++)
      {
         Repart=Repart+real(Help->returnEntry(i, i));
         Impart=Impart+imag(Help->returnEntry(i, i));
      }
      tracer.push_back(Repart);
      tracei.push_back(Impart);
   }
}

#endif

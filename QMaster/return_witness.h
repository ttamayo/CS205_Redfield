#ifndef RETURN_WITNESS_H
#define RETURN_WITNESS_H


#include "headers.h"



class Return_Witness: public ManipulatingObservable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int t_mod;
   bool mum;
   bool left;
   int Ntuples;
   Matrix *mumMatrix;
   vector<double> tracer;
   vector<double> tracei;
private:
   cl_mem d_mum;	// h_ for host d_ for device
   cl_mem d_mup;
   myreal *h_mum;	// help matrix on cpu for mu
   myreal *h_mup;
   int Nsites;
   Matrix *Help;
   Matrix *rho;
   double *h_rhoEntries;
public:
   Return_Witness(System &INsystem, OpenCL_Init &INOpenCLinfo, int INt_mod, bool INmum, bool INleft, int INNtuples, Matrix *INmumMatrix);
   void execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help);
   void freeDeviceMemory();
   ~Return_Witness();
};




Return_Witness::Return_Witness(System &INsystem, OpenCL_Init &INOpenCLinfo, int INt_mod,  bool INmum, bool INleft, int INNtuples, Matrix *INmumMatrix) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo),  t_mod(INt_mod), mum(INmum), left(INleft), Ntuples(INNtuples), mumMatrix(INmumMatrix)
{
   int err;
// copy mum to h_mu
   Nsites=system->Nsites;
   rho=new Matrix(system->Nsites);
   Help=new Matrix(system->Nsites);
   h_rhoEntries=new double[2*system->Nsites*system->Nsites];
   h_mum=new myreal[2*Nsites*Nsites];
   h_mup=new myreal[2*Nsites*Nsites];
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_mum[2*i]=real(system->mum->entry[i]);
      h_mum[2*i+1]=0.;
      h_mup[2*i]=real(system->mup->entry[i]);
      h_mup[2*i+1]=0.;
   }
   d_mum=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*2*Nsites*Nsites, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_mum"<<endl;
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_mum, CL_TRUE, 0, sizeof(myreal)*2*Nsites*Nsites, h_mum,0, NULL, NULL);
   delete [] h_mum;
//
   d_mup=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*2*Nsites*Nsites, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_mup"<<endl;
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_mup, CL_TRUE, 0, sizeof(myreal)*2*Nsites*Nsites, h_mup,0, NULL, NULL);
   delete [] h_mup;
}


Return_Witness::~Return_Witness()
{
   freeDeviceMemory();
}

void Return_Witness::freeDeviceMemory()
{
   clReleaseMemObject(d_mum);
   clReleaseMemObject(d_mup);
}


void Return_Witness::execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help)
{
   if(it%t_mod==0)
   {
      if(mum)		// then I want mum
      {
         Kernel_return_witness((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_sigma_help), (cl_myreal2*)(d_mum), Nsites, left,Ntuples, OpenCLinfo);
      }
      if(!mum)	// then I want mup
      {
         Kernel_return_witness((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_sigma_help), (cl_myreal2*)(d_mup), Nsites, left,Ntuples, OpenCLinfo);
      }
      clEnqueueReadBuffer(OpenCLinfo->queue, d_sigma_help, CL_TRUE, 0, sizeof(double)*2*system->Nsites*system->Nsites, h_rhoEntries,0, NULL, NULL);
      for(int i=0; i<system->Nsites*system->Nsites; i++)
      {
         rho->entry[i]=double_complex(h_rhoEntries[2*i],h_rhoEntries[2*i+1]);
      }
      MatrixMuliplication(Help, mumMatrix, rho, system->Nsites); // first=second * third
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




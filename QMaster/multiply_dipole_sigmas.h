#ifndef MULTIPLY_DIPOLE_SIGMAS_H
#define MULTIPLY_DIPOLE_SIGMAS_H


#include "headers.h"



class Multiply_Dipole_Sigmas: public ManipulatingObservable
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   bool mum;
   bool left;
   int delay;
   int Ntuples;
private:
   cl_mem d_mum;	// h_ for host d_ for device
   cl_mem d_mup;
   myreal *h_mum;	// help matrix on cpu for mu
   myreal *h_mup;
   int Nsites;

public:
   Multiply_Dipole_Sigmas(System &INsystem, OpenCL_Init &INOpenCLinfo, bool INmum, bool INleft, int INdelay, int INNtuples);
   void execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help);
   void freeDeviceMemory();
   ~Multiply_Dipole_Sigmas();
};

Multiply_Dipole_Sigmas::Multiply_Dipole_Sigmas(System &INsystem, OpenCL_Init &INOpenCLinfo, bool INmum, bool INleft, int INdelay, int INNtuples) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), mum(INmum), left(INleft), delay(INdelay), Ntuples(INNtuples)
{
   int err;
// copy mum to h_mu
   Nsites=system->Nsites;
   h_mum=new myreal[2*Nsites*Nsites];
   h_mup=new myreal[2*Nsites*Nsites];
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_mum[2*i]=real(system->mum->entry[i]);
      h_mum[2*i+1]=imag(system->mum->entry[i]);
      h_mup[2*i]=real(system->mup->entry[i]);
      h_mup[2*i+1]=imag(system->mup->entry[i]);
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

Multiply_Dipole_Sigmas::~Multiply_Dipole_Sigmas()
{
   freeDeviceMemory();
}

// provides the different choices of multiplications with dipole matrices!
void Multiply_Dipole_Sigmas::execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help)
{
   if (it==delay)
   {
      if(mum)		// then I want mum
      {
         if(Nsites<=32)
         {
            Kernel_MuMultiplication((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_mum), Nsites, left,Ntuples, OpenCLinfo);		// multiply all sigmas
         }
         else
         {
            Kernel_MuMultiplicationLargeSystem((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_sigma_help), (cl_myreal2*)(d_mum), Nsites, left,Ntuples, OpenCLinfo);
         }
      }
      if(!mum)	// then I want mup
      {
         if(Nsites<=32)
         {
            Kernel_MuMultiplication((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_mup), Nsites, left,Ntuples, OpenCLinfo);		// multiply all sigmas
         }
         else
         {
            Kernel_MuMultiplicationLargeSystem((cl_myreal2*)(d_sigma_System), (cl_myreal2*)(d_sigma_help), (cl_myreal2*)(d_mup), Nsites, left,Ntuples, OpenCLinfo);
         }
      }
   }
}

void Multiply_Dipole_Sigmas::freeDeviceMemory()
{
   clReleaseMemObject(d_mum);
   clReleaseMemObject(d_mup);
}

#endif

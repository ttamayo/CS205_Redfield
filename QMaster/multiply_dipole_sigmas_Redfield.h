#ifndef Multiply_Dipole_Sigmas_Redfield_REDFIELD_H
#define Multiply_Dipole_Sigmas_Redfield_REDFIELD_H


#include "headers.h"



class Multiply_Dipole_Sigmas_Redfield: public ManipulatingObservable
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
   Matrix* h_mum_eig; //mum rotated in eigenbasis
   Matrix* h_mup_eig; // mup rotated in eigenbasis
   myreal *h_mum;	// help matrix on cpu for mu
   myreal *h_mup;
   int Nsites;

public:
   Multiply_Dipole_Sigmas_Redfield(System &INsystem, OpenCL_Init &INOpenCLinfo, bool INmum, bool INleft, int INdelay, int INNtuples);
   void execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help);
   void freeDeviceMemory();
   ~Multiply_Dipole_Sigmas_Redfield();
};

Multiply_Dipole_Sigmas_Redfield::Multiply_Dipole_Sigmas_Redfield(System &INsystem, OpenCL_Init &INOpenCLinfo, bool INmum, bool INleft, int INdelay, int INNtuples) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), mum(INmum), left(INleft), delay(INdelay), Ntuples(INNtuples)
{
   int err;
// copy mum to h_mu
   Nsites=system->Nsites;
   h_mum_eig=new Matrix(Nsites);
   h_mup_eig=new Matrix(Nsites);
   h_mum_eig->assignEntries(*system->mum);
   h_mup_eig->assignEntries(*system->mup);
   system->convertEigenbasis(h_mum_eig);
   system->convertEigenbasis(h_mup_eig);
   h_mum=new myreal[2*Nsites*Nsites];
   h_mup=new myreal[2*Nsites*Nsites];
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_mum[2*i]=real(h_mum_eig->entry[i]);
      h_mum[2*i+1]=imag(h_mum_eig->entry[i]);
      h_mup[2*i]=real(h_mup_eig->entry[i]);
      h_mup[2*i+1]=imag(h_mup_eig->entry[i]);;
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

Multiply_Dipole_Sigmas_Redfield::~Multiply_Dipole_Sigmas_Redfield()
{
   freeDeviceMemory();
}

// provides the different choices of multiplications with dipole matrices!
void Multiply_Dipole_Sigmas_Redfield::execute(double time, int it, cl_mem d_sigma_System, cl_mem d_sigma_help)
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

void Multiply_Dipole_Sigmas_Redfield::freeDeviceMemory()
{
   clReleaseMemObject(d_mum);
   clReleaseMemObject(d_mup);
}

#endif

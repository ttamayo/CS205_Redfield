#ifndef LIOUVILLEHEXCITONREDFIELD_H
#define LIOUVILLEHEXCITONREDFIELD_H


#include "headers.h"


class LiouvilleHExcitonRedfield: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   myreal hbar;
   cl_mem d_Ham;
   myreal *h_Ham;
private:
   int Nsites;
public:
   LiouvilleHExcitonRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo);
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory();
};


LiouvilleHExcitonRedfield::LiouvilleHExcitonRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo)
{
   hbar=system->hbar;
   Nsites=system->Nsites;
   int err;

// Hamiltonian, Redfield propagates in Eigenbasis
   h_Ham=new myreal[Nsites*Nsites*2];
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_Ham[2*i]  =(myreal)real(system->HamDiag->entry[i]);
      h_Ham[2*i+1]=(myreal)imag(system->HamDiag->entry[i]);
   }

   d_Ham=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Nsites*Nsites*2, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Ham"<<endl;
      exit(1);
   }

   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Ham, CL_TRUE, 0, sizeof(myreal)*2*Nsites*Nsites,h_Ham,0, NULL, NULL);
   //  CL_TRUE ensures  that function waits until memcopy finished (equivalent to thread_syncronise), 0, NULL, NULL says that there
   // is no wait list and no event we need to wait for
}

///EVALUATION IN SITEBASIS
void LiouvilleHExcitonRedfield::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   if(Nsites<=32)
   {
      Kernel_HexcitonRedfield(INsigmanew, INsigmaold, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
   else
   {
      Kernel_HexcitonRedfieldLargeSystem(INsigmanew, INsigmaold, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
}
void LiouvilleHExcitonRedfield::freeDeviceMemory()
{
   clReleaseMemObject(d_Ham);
   delete[] h_Ham;
}

#endif

#ifndef LIOUVILLEHEXCITON_H
#define LIOUVILLEHEXCITON_H


#include "headers.h"


class LiouvilleHExciton: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Ntuples; //total numbers of sigma matrices
   myreal hbar;
   cl_mem d_Ham;
   myreal *h_Ham;
private:
   int Nsites;
public:
   LiouvilleHExciton(System &INsystem, OpenCL_Init &INOpenCLinfo, int INNtuples);
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory();
};


LiouvilleHExciton::LiouvilleHExciton(System &INsystem, OpenCL_Init &INOpenCLinfo, int INNtuples) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), Ntuples(INNtuples)
{
   hbar=system->hbar;
   Nsites=system->Nsites;
   int err;

// Hamiltonian
   h_Ham=new myreal[Nsites*Nsites*2];
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_Ham[2*i]  =(myreal)real(system->Ham->entry[i]);
      h_Ham[2*i+1]=(myreal)imag(system->Ham->entry[i]);
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
void LiouvilleHExciton::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   if(Nsites<=32)
   {
      Kernel_Hexciton(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
   else
   {
      Kernel_HexcitonLargeSystem(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
}
void LiouvilleHExciton::freeDeviceMemory()
{
   clReleaseMemObject(d_Ham);
   delete[] h_Ham;
}

#endif

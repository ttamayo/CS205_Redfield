/*
 * Liouville_antiHermitian.h
 *
 *  Created on: Dec 17, 2013
 *      Author: christoph
 */

#ifndef LIOUVILLE_ANTIHERMITIAN_H_
#define LIOUVILLE_ANTIHERMITIAN_H_

#include "headers.h"

class Liouville_antiHermitian: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Ntuples; //total numbers of sigma matrices
   myreal hbar;
   cl_mem d_Gamma; // anti hermitian hamiltonian
   myreal *h_Gamma;
private:
   int Nsites;
public:
   Liouville_antiHermitian(System &INsystem, OpenCL_Init &INOpenCLinfo,  int INNtuples);
   void add_sink(vector<int> site_to_sink, vector<double> rate_to_sink);
   void initialize(void);
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory(void);
};

Liouville_antiHermitian::Liouville_antiHermitian(System &INsystem, OpenCL_Init &INOpenCLinfo, int INNtuples) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), Ntuples(INNtuples)
{
   hbar=system->hbar;
   Nsites=system->Nsites;
// Hamiltonian
   h_Gamma=new myreal[Nsites*Nsites*2];
   for(int i=0; i<Nsites*Nsites*2; i++)
   {
      h_Gamma[i]=0.;
   }
}

void Liouville_antiHermitian::add_sink(vector<int> site_to_sink, vector<double> rate_to_sink)
{
   for(int i=0; i<site_to_sink.size(); i++)
   {
      int id_ii=site_to_sink[i]*Nsites+site_to_sink[i];
      h_Gamma[2*id_ii]+=(myreal)(hbar/rate_to_sink[i]/2);
      h_Gamma[2*id_ii+1]+=(myreal)(0.);
   }
}

void Liouville_antiHermitian::initialize()
{
   int err;
   cl_uint memsize=Nsites*Nsites*2; //factor 2 is needed to get dim2 myreal which represents complex number
   d_Gamma=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Gamma"<<endl;
      exit(1);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Gamma, CL_TRUE, 0, sizeof(myreal)*memsize,h_Gamma,0, NULL, NULL);
}


//EVALUATION IN SITEBASIS
void Liouville_antiHermitian::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   if(Nsites<=32)
   {
      Kernel_antiHermitian(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Gamma), OpenCLinfo);
   }
   else
   {
      Kernel_antiHermitianLargeSystem(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Gamma), OpenCLinfo);
   }
}

void Liouville_antiHermitian::freeDeviceMemory()
{
   clReleaseMemObject(d_Gamma);
   delete[] h_Gamma;
}


#endif /* GPULIOUVILLE_ANTIHERMITIAN_H_ */

#ifndef LIOUVILLEHPULSE_H
#define LIOUVILLEHPULSE_H


#include "headers.h"


class LiouvilleHPulse: public LiouvilleTime
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Ntuples; //total numbers of sigma matrices
   double sigma_pulse; //with of the gaussian pulse in units of seconds
   double center_pulse; //gaussian pulse centered in time around center_pulse in s
   double strength_pulse; //pulse strengt in cm^1/Debye
   double omega_pulse; // carrier frequency of the pulse in cm^{-1}
   myreal hbar;
   cl_mem d_Ham;
   myreal *h_Ham;
private:
   int Nsites;
public:
   LiouvilleHPulse(System &INsystem, OpenCL_Init &INOpenCLinfo, int INNtuples, double INsigma_pulse, double INcenter_pulse, double INstrength_pulse, double INomega_pulse); /// FIXME add Pulse Parameter
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt, myreal time);
   double_complex electricField(double time);
   void freeDeviceMemory();
};

LiouvilleHPulse::LiouvilleHPulse(System &INsystem,  OpenCL_Init &INOpenCLinfo, int INNtuples, double INsigma_pulse, double INcenter_pulse, double INstrength_pulse, double INomega_pulse) :
   system(&INsystem),
   OpenCLinfo(&INOpenCLinfo),
   Ntuples(INNtuples),
   sigma_pulse(INsigma_pulse),
   center_pulse(INcenter_pulse),
   strength_pulse(INstrength_pulse),
   omega_pulse(INomega_pulse)
{
   hbar=system->hbar;
   Nsites=system->Nsites;
   int err;
// Hamiltonian
   h_Ham=new myreal[Nsites*Nsites*2];
   system->create_mu();
// reserve memory for Device
   d_Ham=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*2*Nsites*Nsites, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Ham Pulse"<<endl;
      exit(1);
   }
}

double_complex LiouvilleHPulse::electricField(double time) // pulse strength in cm^1
{
   double_complex amplitude = strength_pulse/ sqrt(2 * M_PI)*exp(double_complex(0,1.)*omega_pulse/hbar*time); // strength_pulse in cm^-1
   return amplitude*exp(-1./(2*sigma_pulse*sigma_pulse)*(time-center_pulse)*(time-center_pulse));
}



void LiouvilleHPulse::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt, myreal time)
{
/// FIXME BEGIN: specify time dependent Hamiltonian
/// INtime is in units of s,  H(t)=E(t)*mu
   Matrix help(Nsites);
   for(int i=1; i<Nsites; i++)
   {
      help.assignEntry(0,i,system->DipoleOperator[i]*electricField(time));
      help.assignEntry(i,0,conj(system->DipoleOperator[i]*electricField(time)));
   }
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_Ham[2*i]=real(help.entry[i]);
      h_Ham[2*i+1]=imag(help.entry[i]);
   }

   ///FIXME how to handel complex parts?
   // copy to Device Memory
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Ham, CL_TRUE, 0, sizeof(myreal)*2*Nsites*Nsites,h_Ham,0, NULL, NULL);
   //  CL_TRUE ensures  that function waits until memcopy finished (equivalent to thread_syncronise), 0, NULL, NULL says that there
   // Propagate time dependent Hamiltonian in the same way as the exciton Hamiltonian
   if(Nsites<=32)
   {
      Kernel_Hexciton(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
   else
   {
      Kernel_HexcitonLargeSystem(INsigmanew, INsigmaold, Ntuples, Nsites, hbar,dt, (cl_myreal2*)(d_Ham), OpenCLinfo);
   }
}

void LiouvilleHPulse::freeDeviceMemory()
{
   clReleaseMemObject(d_Ham);
   delete[] h_Ham;
}



#endif

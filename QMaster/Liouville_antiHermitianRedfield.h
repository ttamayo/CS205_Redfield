/*
 * Liouville_antiHermitianRedfield.h
 *
 *  Created on: May 03, 2014
 *      Author: christoph
 */

#ifndef LIOUVILLE_ANTIHERMITIANREDFIELD_H_
#define LIOUVILLE_ANTIHERMITIANREDFIELD_H_

#include "headers.h"

class Liouville_antiHermitianRedfield: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   int Ntuples; //total numbers of sigma matrices
   myreal hbar;
   myreal *h_Gamma;
   Matrix *rho_new;
   Matrix *rho_old;
   Matrix *Gamma;
   myreal *h_rhoEntries_new;
   myreal *h_rhoEntries_old;
private:
   int Nsites;
public:
   Liouville_antiHermitianRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo,  int INNtuples);
   void add_sink(vector<int> site_to_sink, vector<double> rate_to_sink);
   void initialize(void);
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
   void freeDeviceMemory(void);
};

Liouville_antiHermitianRedfield::Liouville_antiHermitianRedfield(System &INsystem, OpenCL_Init &INOpenCLinfo, int INNtuples) :
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

   rho_new=new Matrix(system->Nsites);
   rho_old=new Matrix(system->Nsites);
   h_rhoEntries_new=new double[2*system->Nsites*system->Nsites];
   h_rhoEntries_old=new double[2*system->Nsites*system->Nsites];

}

void Liouville_antiHermitianRedfield::add_sink(vector<int> site_to_sink, vector<double> rate_to_sink)
{
   for(int i=0; i<site_to_sink.size(); i++)
   {
      int id_ii=site_to_sink[i]*Nsites+site_to_sink[i];
      //multiply with im
      h_Gamma[2*id_ii]+=(myreal)(hbar/rate_to_sink[i]/2);
      h_Gamma[2*id_ii+1]+=(myreal)(0.);
   }
}

void Liouville_antiHermitianRedfield::initialize()
{
   Gamma=new Matrix(Nsites);
   for(int i=0; i<Nsites*Nsites; i++)
   {
      Gamma->entry[i]=double_complex(h_Gamma[2*i],h_Gamma[2*i+1]);
   }
}


//EVALUATION IN SITEBASIS
void Liouville_antiHermitianRedfield::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
	 clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);
	 clEnqueueReadBuffer(OpenCLinfo->queue, (cl_mem) (INsigmaold), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_old,0, NULL, NULL);


	 for(int i=0; i<system->Nsites*system->Nsites; i++)
	    {
	       rho_new->entry[i]=double_complex(h_rhoEntries_new[2*i],h_rhoEntries_new[2*i+1]);
	       rho_old->entry[i]=double_complex(h_rhoEntries_old[2*i],h_rhoEntries_old[2*i+1]);
	    }
    //do anti-hermitian in site-basis (does not work in eigenbasis!
	 system->convertSitebasis(rho_old);
	 Matrix Help1(Nsites);
	 // Gamma*rho
	 MatrixMuliplication(&Help1,Gamma, rho_old, Nsites);
	 Matrix Help2(Nsites);
	 // rho*Gamma
	 MatrixMuliplication(&Help2, rho_old, Gamma, Nsites);
	 //Gamma*rho+rho*Gamma
	 MatrixAdd(&Help1, &Help2, Nsites);
	 //*-dt/hbar
	 MultiplyScalar(&Help1, -dt/system->hbar, Nsites);
	 //
	 //go back to eigenbasis
	 system->convertEigenbasis(&Help1);
	 MatrixAdd(rho_new, &Help1, Nsites);
	   for(int i=0; i<system->Nsites*system->Nsites; i++)
	   {
	      h_rhoEntries_new[2*i]=real(rho_new->entry[i]);
	      h_rhoEntries_new[2*i+1]=imag(rho_new->entry[i]);
	   }
	 clEnqueueWriteBuffer(OpenCLinfo->queue, (cl_mem) (INsigmanew), CL_TRUE, 0, sizeof(myreal)*2*system->Nsites*system->Nsites,h_rhoEntries_new,0, NULL, NULL);
}

void Liouville_antiHermitianRedfield::freeDeviceMemory()
{
   delete[] h_Gamma;
}


#endif /* GPULIOUVILLE_ANTIHERMITIANREDFIELD_H_ */

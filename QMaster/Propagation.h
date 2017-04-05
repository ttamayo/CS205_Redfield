#ifndef PROPAGATION
#define PROPAGATION

#include "headers.h"

class Propagation
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   SigmaHelpMatrices *sigma_Host;
   int Nsteps;
   double dt;
   int Nmax;
   int Nsites;
   int Ncoupled;
   int Ntuples;
   myreal *h_rhoEntries;
   bool breakProp;
   // Sigma Matrices, real System
   cl_mem d_sigma_System;
private:
   cl_int err;
   double time;
   vector<Liouville*> All_Liouville;
   vector<LiouvilleTime*> All_LiouvilleTime;
   vector<Observable*> All_Observables;
   vector<ManipulatingObservable*> All_ManipulatingObservables;
   vector<BreakRule*>  All_BreakRules; //defines rules for breaking the propagation
   vector<double> breakCond;  //condition for Breaking the propagation
   double* breakVal;  // actual value, if this reaches breakCond then propagation stops
   cl_mem d_sigma_k1; //FIXME versuche mit weniger sigma_k's auszukommen
   cl_mem d_sigma_k2;
   cl_mem d_sigma_k3;
   cl_mem d_sigma_k4;
   cl_mem d_sigma_help;
public:
   Propagation(System &INsystem, OpenCL_Init &INOpenCLinfo, SigmaHelpMatrices &INsigma_Host, int INNsteps, double INdt, int INNmax);
   void Add_Liouville(Liouville &INLiouville);
   void Add_LiouvilleTime(LiouvilleTime &INGPULiouvilleTime);
   void Runge_Kutta(cl_myreal2* INsigma_System, double dt, double time);
   void propagate(int Nmod);
   void initialize(Matrix rhoinit);
   void Add_Observable(Observable &INObservable);
   void Add_ManipulatingObservable(ManipulatingObservable &INManipulatingObservable);
   void freeDeviceMemory();
   void clear_Sigma();
   void clear_ManipulatingObservable();
   void Add_BreakRule(BreakRule &INBreakRule, double breakcond, int id);
};


Propagation::Propagation(System &INsystem,  OpenCL_Init &INOpenCLinfo, SigmaHelpMatrices &INsigma_Host, int INNsteps, double INdt, int INNmax) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), sigma_Host(&INsigma_Host), Nsteps(INNsteps), dt(INdt), Nmax(INNmax)
{
   breakProp=false; //initialize bool for break conditions for time evolution
   breakVal=new double[10]; //allow maximal 10 break conditions
   Nsites=system->Nsites;
   Ncoupled=system->Ncoupled;
   h_rhoEntries=new myreal[2*Nsites*Nsites];
   Ntuples=sigma_Host->AllsigmaTuples.size();
/// Allocated GPU Memory
//d_Memory of all used sigmamatrices
   unsigned long int memsize=(unsigned long int) (Nsites*Ntuples*2)*(unsigned long int)(Nsites); //factor 2 is needed to get dim2 myreal which represents complex number
   fprintf(stdout,"# sigma matrices, start allocate DEVICE-Memory: %f MB \n", 6.*sizeof(myreal)*memsize/1024./1024.);
   d_sigma_k1=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_k1 (out of memory)\n");
      exit(0);
   }
   d_sigma_k2=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_k2 (out of memory)\n");
      exit(0);
   }
   d_sigma_k3=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_k3 (out of memory)\n");
      exit(0);
   }
   d_sigma_k4=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_k4 (out of memory)\n");
      exit(0);
   }
   d_sigma_help=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_help (out of memory)\n");
      exit(0);
   }
   d_sigma_System=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_WRITE, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      fprintf(stderr,"Error: Couldn't allocate d_sigma_System (out of memory)\n");
      exit(0);
   }
}

void Propagation::clear_Sigma()
{
   InitializeSigma((cl_myreal2*) d_sigma_k1, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k2, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k3, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k4, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_help, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_System, Ntuples, Nsites,OpenCLinfo);
}

void Propagation::clear_ManipulatingObservable()
{
   All_ManipulatingObservables.erase(All_ManipulatingObservables.begin(), All_ManipulatingObservables.end() );
}
void Propagation::initialize(Matrix rhoinit)
{
   InitializeSigma((cl_myreal2*) d_sigma_System, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*)  d_sigma_help, Ntuples, Nsites, OpenCLinfo);
//initialize rho_0..0 on GPU with h_rho
   for(int i=0; i<Nsites*Nsites; i++)
   {
      h_rhoEntries[2*i]=real(rhoinit.entry[i]);
      h_rhoEntries[2*i+1]=imag(rhoinit.entry[i]);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_sigma_System, CL_TRUE, 0, sizeof(myreal)*2*Nsites*Nsites,h_rhoEntries,0, NULL, NULL);    //  CL_TRUE ensures  that function waits until memcopy finished (equivalent to thread_syncronise), 0, NULL, NULL says that there is no wait list and no event we need to wait for
}

void Propagation::Add_Liouville(Liouville &INLiouville)
{
   All_Liouville.push_back(&INLiouville);
}

void Propagation::Add_LiouvilleTime(LiouvilleTime &INLiouvilleTime)
{
   All_LiouvilleTime.push_back(&INLiouvilleTime);
}

void Propagation::Add_ManipulatingObservable(ManipulatingObservable &INManipulatingObservable)
{
   All_ManipulatingObservables.push_back(&INManipulatingObservable);
}

void Propagation::Add_Observable(Observable &INObservable)
{
   All_Observables.push_back(&INObservable);
}

void Propagation::Add_BreakRule(BreakRule &INBreakRule, double breakcond, int id)
{
   All_BreakRules.push_back(&INBreakRule);
   breakCond.push_back(breakcond);
   //initialize breakval
   breakVal[id]=0.;
}

void Propagation::Runge_Kutta(cl_myreal2* INsigma_System, double dt, double time)
{
//Reset d_sigma_k1, d_sigma_k2, d_sigma_k3, d_sigma_k4
   InitializeSigma((cl_myreal2*) d_sigma_k1, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k2, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k3, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_k4, Ntuples, Nsites,OpenCLinfo);
   InitializeSigma((cl_myreal2*) d_sigma_help, Ntuples, Nsites,OpenCLinfo);

/// compute d_sigma_k1 = dt *Liouv(INsigma_System)
   for(int j=0; j<All_Liouville.size(); j++)
   {
      (*All_Liouville[j]).execute((cl_myreal2*) INsigma_System, (cl_myreal2*) d_sigma_k1, dt); //Add contribution of GPULiouville[j]
   }
   for(int j=0; j<All_LiouvilleTime.size(); j++)
   {
      (*All_LiouvilleTime[j]).execute((cl_myreal2*) INsigma_System, (cl_myreal2*) d_sigma_k1, dt, time); //Add contribution of GPULiouville[j]
   }
//
/// compute d_sigma_help=1/2*d_sigma_k1+INsigma_system
   aXsig1Plussig2((cl_myreal2*)d_sigma_help, (cl_myreal2*) d_sigma_k1, INsigma_System, 0.5f, Ntuples, Nsites,OpenCLinfo);
//
/// compute d_sigma_k2 = dt *Liouv(d_sigma_help)
   for(int j=0; j<All_Liouville.size(); j++)
   {
      (*All_Liouville[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k2, dt); //Add contribution of GPULiouville[j]
   }
   for(int j=0; j<All_LiouvilleTime.size(); j++)
   {
      (*All_LiouvilleTime[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k2, dt, time+0.5*dt); //Add contribution of GPULiouville[j]
   }
//
/// compute d_sigma_help=1/2*d_sigma_k2+INsigma_system
   aXsig1Plussig2((cl_myreal2*)d_sigma_help, (cl_myreal2*) d_sigma_k2, INsigma_System, 0.5f, Ntuples, Nsites,OpenCLinfo);
//
/// compute d_sigma_k3 = dt *Liouv(d_sigma_help)
   for(int j=0; j<All_Liouville.size(); j++)
   {
      (*All_Liouville[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k3, dt); //Add contribution of GPULiouville[j]
   }
   for(int j=0; j<All_LiouvilleTime.size(); j++)
   {
      (*All_LiouvilleTime[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k3, dt, time+0.5*dt); //Add contribution of GPULiouville[j]
   }
//
/// compute d_sigma_help=d_sigma_k3+INsigma_system
   aXsig1Plussig2((cl_myreal2*)d_sigma_help, (cl_myreal2*) d_sigma_k3, INsigma_System, 1.f, Ntuples, Nsites,OpenCLinfo);
//
/// compute d_sigma_k4 = dt *Liouv(d_sigma_help)
   for(int j=0; j<All_Liouville.size(); j++)
   {
      (*All_Liouville[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k4, dt); //Add contribution of GPULiouville[j]
   }
   for(int j=0; j<All_LiouvilleTime.size(); j++)
   {
      (*All_LiouvilleTime[j]).execute((cl_myreal2*) d_sigma_help, (cl_myreal2*) d_sigma_k4, dt, time+dt); //Add contribution of GPULiouville[j]
   }
//
//
/// compute INsigma_System=INsigma_System+1/6*d_sigma_k1+1/3*d_sigma_k2+1/3*d_sigma_k3+1/6*d_sigma_k4
   ReturnNewSigma(INsigma_System, (cl_myreal2*) d_sigma_k1, (cl_myreal2*) d_sigma_k2, (cl_myreal2*) d_sigma_k3, (cl_myreal2*) d_sigma_k4, Ntuples, Nsites,OpenCLinfo);
}

void Propagation::propagate(int Nmod)
{
   time=0.;
// evaluate all observables
   for(int j=0; j<All_Observables.size(); j++)
   {
      (*All_Observables[j]).execute(time, 0, d_sigma_System);
   }
   for(int j=0; j<All_ManipulatingObservables.size(); j++)
   {
      (*All_ManipulatingObservables[j]).execute(time, 0, d_sigma_System, d_sigma_help);
   }
   //evaluate break rule
   for(int j=0; j<All_BreakRules.size(); j++)
   {
      (*All_BreakRules[j]).execute(breakVal, d_sigma_System,0);
   }
///time evolution
   for(int k=0; k<Nsteps; k++)
   {
      time=(k+1)*dt; // time for new rho(t), for time dependent hamiltonian need to evaluate at previous timestep t-dt
      if((k+1)%Nmod==0)
      {
         //  progress meter
         cout<<"# PROGRESS=>"<<100*(k+1)/Nsteps<<" percent done"<<endl;
         cout.flush();
      }
      Runge_Kutta((cl_myreal2*) d_sigma_System, dt, time-dt);// time for new rho(t), for time dependent hamiltonian need to evaluate at previous timestep t-dt
      // evaluate all observables
      for(int j=0; j<All_Observables.size(); j++)
      {
         (*All_Observables[j]).execute(time, k+1, d_sigma_System);
      }
      for(int j=0; j<All_ManipulatingObservables.size(); j++)
      {
         (*All_ManipulatingObservables[j]).execute(time, k+1, d_sigma_System,d_sigma_help);
      }
      // Check break rules
      for(int j=0; j<All_BreakRules.size(); j++)
      {
         (*All_BreakRules[j]).execute(breakVal, d_sigma_System, k+1);
         if(breakVal[j]<breakCond[j])
         {
            breakProp=true;
         }
      }
      if(breakProp==true)
      {
         cout<<"# reached one of the break conditions:"<<endl;
         for(int j=0; j<All_BreakRules.size(); j++)
         {
            cout<<"# breakCond "<<j<<" :  breakVal="<<scientific<<breakVal[j]<<"  breakCond="<<scientific<<breakCond[j]<<endl;
         }
         cout << scientific;
         cout<<"# break at it="<<k+1<<endl;
         break;
      }
   }
   if(All_BreakRules.size()>0 and breakProp==false)
   {
      cout<<"# warning: break conditions are not yet fulfilled:"<<endl;
      cout<<"# warning: reached end of Nsteps="<<Nsteps<<"; break anyway!"<<endl;
   }
}


void Propagation::freeDeviceMemory()
{
   cout<<"# free sigma matrices"<<endl;
   clReleaseMemObject(d_sigma_k1);
   clReleaseMemObject(d_sigma_k2);
   clReleaseMemObject(d_sigma_k3);
   clReleaseMemObject(d_sigma_k4);
   clReleaseMemObject(d_sigma_help);
   clReleaseMemObject(d_sigma_System);
}

#endif

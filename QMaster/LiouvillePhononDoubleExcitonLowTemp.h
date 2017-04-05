#ifndef LIOUVILLEPHONONDoubleExcitonLowTemp_H
#define LIOUVILLEPHONONDoubleExcitonLowTemp_H

#include "headers.h"

class LiouvillePhononDoubleExcitonLowTemp: public Liouville
{
public:
   System *system;
   OpenCL_Init* OpenCLinfo;
   SigmaHelpMatrices *sigma_Host;
   int Nmax;
   double hbar;
   int Nsites;
   int Nsites1;
   int Ncoupled; //number of sites coupled to phonons
   int NsitesCoupled;
private:
   vector<int> Sites;
   vector<double> listlambda;
   vector<double_complex> listgamma;
   vector<int> listNDL;
   vector<int> listNDLsite;
   cl_mem d_TupleNorm;
   int* h_Sites;
   myreal *h_listlambda;
   myreal *h_listgamma;
   cl_mem d_Sites;
   cl_mem d_listlambda;
   cl_mem d_listgamma;
   int* h_listNDL; //Number of Drude-Lorentz peaks per sites
   int* h_listNDLsite;//actual site number
   cl_mem d_listNDL; //Number of Drude-Lorentz peaks per sites
   cl_mem d_listNDLsite;//actual site number
   cl_mem d_AllSigmaTuples; //Defines Memory for all sigma Tuples, maximal Order Nmax
   myreal *h_Tupleprefact; //prefact is defines diagonal coupling term n_11*gamma_11+n_12*gamma_12 ...
   cl_mem d_Tupleprefact;
   cl_mem d_MemIDnplus1; //Speicheraddresse der Sigmamatrizen zu n+1
   int *h_Memnplus1start;
   cl_mem d_Memnplus1start;
   cl_mem d_MemIDnminus1; //Speicheraddresse der Sigmamatrizen zu n-1
   myreal *h_prefactLowTemp;
   cl_mem d_prefactLowTemp;
   cl_mem d_prefactLowTemp2;
   myreal *h_prefactLowTemp2;
   myreal *h_S0;
   cl_mem d_S0;
   myreal *h_Imchi0;
   cl_mem d_Imchi0;
   myreal *h_Rechi0;
   cl_mem d_Rechi0;
   int Ntuples;
public:
   LiouvillePhononDoubleExcitonLowTemp(System &INsystem, OpenCL_Init &INOpenCLinfo, SigmaHelpMatrices &INsigma_Host, int INNmax);
   void AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma); //INNDL stands for the total number of Drude-Lorentz peaks in the spectral density for site INsite
   void initialize();
   void freeDeviceMemory();
   void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt);
};

LiouvillePhononDoubleExcitonLowTemp::LiouvillePhononDoubleExcitonLowTemp(System &INsystem, OpenCL_Init &INOpenCLinfo, SigmaHelpMatrices &INsigma_Host, int INNmax) :
   system(&INsystem), OpenCLinfo(&INOpenCLinfo), sigma_Host(&INsigma_Host), Nmax(INNmax)
{
   hbar=system->hbar;
   Nsites=system->Nsites;
   Nsites1=system->Nsites1;
   Ntuples=sigma_Host->AllsigmaTuples.size();
}

void LiouvillePhononDoubleExcitonLowTemp::AddSite(int INsite, int INNDL, double *INlambda, double_complex *INgamma)
{
   listNDL.push_back(INNDL);
   listNDLsite.push_back(INsite);
   for(int i=0; i<INNDL; i++)
   {
      Sites.push_back(INsite);
      listlambda.push_back(INlambda[i]);
      listgamma.push_back(INgamma[i]);
//   cout<<"site "<<INsite<<" lambda "<<INlambda[i]/system->invcmtomeV<<" invgamma in fs "<<1.e15/INgamma[i]<<endl;
   }
}

void LiouvillePhononDoubleExcitonLowTemp::initialize()
{
   unsigned long int totalmemsize=0;
   Ncoupled=Sites.size();
   int err;
   NsitesCoupled=listNDLsite.size();
   if(Ncoupled!=system->Ncoupled)
   {
      fprintf(stderr,"###\n");
      fprintf(stderr,"###\n");
      fprintf(stderr,"### ERROR: NUMBER OF COUPLED SITES IN GPULIOUVILLEPHONON LINE 76 DOES NOT COINCIDE WITH NCOUPLED DEFINED IN SYSTEM.H\n");
      fprintf(stderr,"### GPULiouvillePhonon.h: Ncoupled= %i and System.h: Ncoupled= %i\nABORT.\n", Ncoupled, system->Ncoupled);
      exit(1);
   }
   h_Sites=new int[Ncoupled];
   h_listgamma=new double[2*Ncoupled]; //gamma complex, need memory for real and imaginary part
   h_listlambda=new double[Ncoupled];
   h_prefactLowTemp=new double[2*Ncoupled];
   double v1;
   v1=2.*M_PI/(system->beta*system->hbar);
   for(int i=0; i<Ncoupled; i++)
   {
      h_Sites[i]=Sites[i];
      h_listgamma[2*i]=real(listgamma[i]);
      h_listgamma[2*i+1]=imag(listgamma[i]);
      h_listlambda[i]=listlambda[i];
      double_complex hilf;
      hilf=-double_complex(0.,1.)*2.*listlambda[i]/(system->beta*system->hbar);
      hilf*=2.*listgamma[i]/(v1*v1-listgamma[i]*listgamma[i])*listgamma[i];
      h_prefactLowTemp[2*i]=real(hilf);
      h_prefactLowTemp[2*i+1]=imag(hilf);
   }
   h_prefactLowTemp2=new double[2*NsitesCoupled];
   int id=0;
   double_complex hilf2=0.;
   for(int j=0; j<NsitesCoupled; j++)
   {
      for(int i=0; i<listNDL[j]; i++)
      {
         hilf2+=4.*listlambda[id]/(system->beta*system->hbar*system->hbar)*listgamma[id]/(v1*v1-listgamma[id]*listgamma[id]);
         id+=1;
      }
      h_prefactLowTemp2[2*j]=-real(hilf2);
      h_prefactLowTemp2[2*j+1]=-imag(hilf2);
      hilf2=0.;
   }
   // GPU-init: d_Sites
   totalmemsize+=sizeof(int)*Ncoupled;
   totalmemsize+=sizeof(int)*NsitesCoupled;
   totalmemsize+=sizeof(int)*NsitesCoupled;
   totalmemsize+=sizeof(int)*NsitesCoupled;
   totalmemsize+=sizeof(int)*Ntuples*Ncoupled;
   totalmemsize+=sizeof(int)*Ntuples;
   totalmemsize+=sizeof(int)*Ntuples*Ncoupled;
   totalmemsize+=sizeof(int)*Ncoupled*Ntuples;
   totalmemsize+=sizeof(myreal)*Ncoupled*2;
   totalmemsize+=sizeof(myreal)*Ncoupled;
   totalmemsize+=sizeof(myreal)*Ntuples*2;
   totalmemsize+=sizeof(myreal)*Ncoupled;
   totalmemsize+=sizeof(myreal)*Ncoupled;
   totalmemsize+=sizeof(myreal)*Ncoupled;
   totalmemsize+=sizeof(myreal)*Ncoupled*2;
   totalmemsize+=sizeof(myreal)*NsitesCoupled*2;
//
   fprintf(stdout,"# create linked list, allocated CPU-Memory, will be freed again while allocating DEVICE-Memory: %f MB\n", totalmemsize/1024./1024.);
   fprintf(stdout,"# create linked list, allocated DEVICE-Memory: %f MB\n", totalmemsize/1024./1024.);
//
//
//
// GPU-init: d_Sites
   d_Sites=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*Ncoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Sites (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Sites, CL_TRUE, 0, sizeof(int)*Ncoupled, h_Sites,0, NULL, NULL);
   delete [] h_Sites;
// GPU-init: d_listNDL
   h_listNDL=new int[NsitesCoupled];
   for(int i=0; i<NsitesCoupled; i++)
   {
      h_listNDL[i]=listNDL[i];
   }
   d_listNDL=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*NsitesCoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_listNDL (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_listNDL, CL_TRUE, 0, sizeof(int)*NsitesCoupled, h_listNDL, 0, NULL, NULL);
   // delete [] h_listNDL; needed later => delete later
   // GPU-init: d_listNDLsite
   h_listNDLsite=new int[NsitesCoupled];
   for(int i=0; i<NsitesCoupled; i++)
   {
      h_listNDLsite[i]=listNDLsite[i];
   }
   d_listNDLsite=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*NsitesCoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_listNDLsite (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_listNDLsite, CL_TRUE, 0, sizeof(int)*NsitesCoupled, h_listNDLsite,0, NULL, NULL);
   delete [] h_listNDLsite;
// GPU-init: d_listgamma
// gamma complex -> need 2*Ncoupled memory
   d_listgamma=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled*2, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_listgamma (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_listgamma, CL_TRUE, 0, sizeof(myreal)*Ncoupled*2, h_listgamma,0, NULL, NULL);
   //delete [] h_listgamma; //needed later => delete later
// GPU-init: d_listlambda
   d_listlambda=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_listlambda (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_listlambda, CL_TRUE, 0, sizeof(myreal)*Ncoupled, h_listlambda,0, NULL, NULL);
//   delete [] h_listlambda; //needed later => delete later
// init SigmaTuples
   // init SigmaTuples
   sigma_Host->initListAllSigma();
   d_AllSigmaTuples=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*Ncoupled*Ntuples, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_AllSigmaTuples (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_AllSigmaTuples, CL_TRUE, 0, sizeof(int)*Ncoupled*Ntuples, sigma_Host->ListAllsigmaTuples,0, NULL, NULL);
   delete [] sigma_Host->ListAllsigmaTuples;
//
//d_TupleNorm;
   int memsize2=Ntuples;
   d_TupleNorm=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*memsize2, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_TupleNorm (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_TupleNorm, CL_TRUE, 0, sizeof(int)*memsize2, sigma_Host->TupleNorm, 0, NULL, NULL);
   delete [] sigma_Host->TupleNorm;
   //
//
// define prefact for evaluation of diagonl coupling term in hierarchy
   h_Tupleprefact=new myreal[2*Ntuples]; //h_Tupleprefact is list of complex numbers
   for(int ituple=0; ituple<Ntuples; ituple++)
   {
      double Repart_prefact=0.;
      double Impart_prefact=0.;
      for(int j=0; j<Ncoupled; j++)
      {
         Repart_prefact+=sigma_Host->AllsigmaTuples[ituple]->entry[j]*h_listgamma[2*j];
         Impart_prefact+=sigma_Host->AllsigmaTuples[ituple]->entry[j]*h_listgamma[2*j+1];
      }
      //   cout<<prefact<<endl;
      h_Tupleprefact[2*ituple]=Repart_prefact;
      h_Tupleprefact[2*ituple+1]=Impart_prefact;
   }
//  d_Tupleprefact;
   int memsize=Ntuples*2; //since d_Tupleprefact is complex
   d_Tupleprefact=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*memsize, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Tupleprefact (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Tupleprefact, CL_TRUE, 0, sizeof(myreal)*memsize, h_Tupleprefact,0, NULL, NULL);
   delete [] h_Tupleprefact;
//  d_MemIDnplus1; //Speicheraddresse der Sigmamatrizen zu n+1
   memsize=Ntuples*Ncoupled;
   d_MemIDnplus1=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*memsize, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_MemIDnplus1 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_MemIDnplus1, CL_TRUE, 0, sizeof(int)*memsize,sigma_Host->MemIDnplus1,0, NULL, NULL);
   delete [] sigma_Host->MemIDnplus1;
//
   h_Memnplus1start=new int[NsitesCoupled];
   h_Memnplus1start[0]=0;
   int Ntot=0;
   for(int i=1; i<NsitesCoupled; i++)
   {
      Ntot+=h_listNDL[i-1];
      h_Memnplus1start[i]=Ntot;
   }
   delete [] h_listNDL;
   d_Memnplus1start=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*NsitesCoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Memnplus1start (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Memnplus1start, CL_TRUE, 0, sizeof(int)*NsitesCoupled, h_Memnplus1start,0, NULL, NULL);
   delete [] h_Memnplus1start;
//d_MemIDnminus1; //Speicheraddresse der Sigmamatrizen zu n-1
   memsize=Ntuples*Ncoupled;
   d_MemIDnminus1=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(int)*memsize, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_MemIDnminus1 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_MemIDnminus1, CL_TRUE, 0, sizeof(int)*memsize,sigma_Host->MemIDnminus1,0, NULL, NULL);
   delete [] sigma_Host->MemIDnminus1;
//
// Define Correlations
   h_S0=new double[Ncoupled];
   for(int i=0; i<Ncoupled; i++)
   {
      h_S0[i]=2.*h_listlambda[i]/(system->beta*hbar);
   }
   d_S0=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_S0 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_S0, CL_TRUE, 0, sizeof(myreal)*Ncoupled, h_S0, 0, NULL, NULL);
   delete [] h_S0;

//   system->chi0 only imaginary part is non-zero:
// sufficient to consider imag part of chi0
   h_Imchi0=new double[Ncoupled];
   h_Rechi0=new double[Ncoupled];
   for(int i=0; i<Ncoupled; i++)
   {
      h_Imchi0[i]=-h_listlambda[i]*h_listgamma[2*i]; //imaginary part of chi0 goes with realpart of gamma, chi0=-i*lambda*gamma
      h_Rechi0[i]=h_listlambda[i]*h_listgamma[2*i+1]; //Realpart of chi0 goes with imaginary part of gamm
   }
   delete [] h_listgamma;
   delete [] h_listlambda;
   d_Imchi0=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Imchi0 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Imchi0, CL_TRUE, 0, sizeof(myreal)*Ncoupled, h_Imchi0,0, NULL, NULL);
   delete [] h_Imchi0;
   d_Rechi0=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_Rechi0 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_Rechi0, CL_TRUE, 0, sizeof(myreal)*Ncoupled, h_Rechi0,0, NULL, NULL);
   delete [] h_Rechi0;
   d_prefactLowTemp=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*Ncoupled*2, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_prefactLowTemp (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_prefactLowTemp, CL_TRUE, 0, sizeof(myreal)*Ncoupled*2, h_prefactLowTemp,0, NULL, NULL);
   delete [] h_prefactLowTemp;
   d_prefactLowTemp2=clCreateBuffer(OpenCLinfo->context, CL_MEM_READ_ONLY, sizeof(myreal)*NsitesCoupled*2, 0, &err);
   if(err<0)
   {
      cout<<"Error: Couldn't allocate d_prefactLowTemp2 (out of memory)"<<endl;
      exit(0);
   }
   clEnqueueWriteBuffer(OpenCLinfo->queue, d_prefactLowTemp2, CL_TRUE, 0, sizeof(myreal)*NsitesCoupled*2, h_prefactLowTemp2,0, NULL, NULL);
   delete [] h_prefactLowTemp2;
}


void LiouvillePhononDoubleExcitonLowTemp::execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt)
{
   if(Nsites<=32)
   {
      Kernel_HPhononDoubleLowTemp_Part1(INsigmanew, INsigmaold, (int*) d_AllSigmaTuples, (int*) d_MemIDnminus1, (int*) d_MemIDnplus1, (int*) d_Memnplus1start, (int*) d_TupleNorm, Nsites, NsitesCoupled, Ncoupled, Ntuples, Nmax, dt, hbar, (myreal*) d_S0, (myreal*) d_Imchi0, (myreal*) d_Rechi0, (cl_myreal2*) d_prefactLowTemp, (cl_myreal2*) d_prefactLowTemp2, (int*) d_listNDL, (int*) d_listNDLsite, OpenCLinfo);
      Kernel_HPhononDouble_DoubleExcitationLowTemp(INsigmanew, INsigmaold, (int*) d_AllSigmaTuples, (int*) d_MemIDnminus1, (int*) d_MemIDnplus1, (int*) d_Memnplus1start, (int*) d_TupleNorm, Nsites, Nsites1,  NsitesCoupled,  Ncoupled, Ntuples, Nmax, dt, hbar, (myreal*) d_S0, (myreal*) d_Imchi0, (myreal*) d_Rechi0,  (cl_myreal2*) d_prefactLowTemp, (cl_myreal2*) d_prefactLowTemp2, (int*) d_listNDL, (int*) d_listNDLsite, OpenCLinfo);
      Kernel_HPhononDoubleLowTemp_Part2(INsigmanew, INsigmaold, Ncoupled, Nsites, Ntuples, dt, (cl_myreal2*) (d_Tupleprefact), OpenCLinfo);
   }
   else
   {
      Kernel_HPhononDoubleLowTemp_Part1LargeSystem(INsigmanew, INsigmaold, (int*) d_AllSigmaTuples, (int*) d_MemIDnminus1, (int*) d_MemIDnplus1, (int*) d_Memnplus1start, (int*) d_TupleNorm, Nsites, NsitesCoupled, Ncoupled, Ntuples, Nmax, dt, hbar, (myreal*) d_S0, (myreal*) d_Imchi0, (myreal*) d_Rechi0, (cl_myreal2*) d_prefactLowTemp, (cl_myreal2*) d_prefactLowTemp2, (int*) d_listNDL, (int*) d_listNDLsite, OpenCLinfo);
      Kernel_HPhononDouble_DoubleExcitationLowTempLargeSystem(INsigmanew, INsigmaold, (int*) d_AllSigmaTuples, (int*) d_MemIDnminus1, (int*) d_MemIDnplus1, (int*) d_Memnplus1start, (int*) d_TupleNorm, Nsites, Nsites1,  NsitesCoupled,  Ncoupled, Ntuples, Nmax, dt, hbar, (myreal*) d_S0, (myreal*) d_Imchi0, (myreal*) d_Rechi0,  (cl_myreal2*) d_prefactLowTemp, (cl_myreal2*) d_prefactLowTemp2, (int*) d_listNDL, (int*) d_listNDLsite, OpenCLinfo);
      Kernel_HPhononDoubleLowTemp_Part2LargeSystem(INsigmanew, INsigmaold, Ncoupled, Nsites, Ntuples, dt, (cl_myreal2*) (d_Tupleprefact), OpenCLinfo);
   }
}

void LiouvillePhononDoubleExcitonLowTemp::freeDeviceMemory()
{
   cout<<"# free linked list"<<endl;
   clReleaseMemObject(d_TupleNorm);
   clReleaseMemObject(d_Sites);
   clReleaseMemObject(d_listlambda);
   clReleaseMemObject(d_listgamma);
   clReleaseMemObject(d_listNDL); //Number of Drude-Lorentz peaks per sites
   clReleaseMemObject(d_listNDLsite);//actual site number
   clReleaseMemObject(d_AllSigmaTuples); //Defines Memory for all sigma Tuples, maximal Order Nmax
   clReleaseMemObject(d_Tupleprefact);
   clReleaseMemObject(d_MemIDnplus1); //Speicheraddresse der Sigmamatrizen zu n+1
   clReleaseMemObject(d_Memnplus1start);
   clReleaseMemObject(d_MemIDnminus1); //Speicheraddresse der Sigmamatrizen zu n-1
   clReleaseMemObject(d_S0);
   clReleaseMemObject(d_Imchi0);
   clReleaseMemObject(d_Rechi0);
   clReleaseMemObject(d_prefactLowTemp);
   clReleaseMemObject(d_prefactLowTemp2);
}

#endif






#ifndef KERNELSOURCE_H
#define KERNELSOURCE_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;

class Kernel_Source
{
public:
   int numFiles;
   const char** program_buffer;
   size_t* program_size;
public:
   Kernel_Source(void);
   void define_InitSigma(void);
   void define_RungeKutta(void);
   void define_Hexciton(void);
   void define_H_antiHermitian(void);
   void define_HPhononLowTemp(void);
   void define_HPhononDoubleExcitonLowTemp(void);
   void define_multiply_dipole_sigmas(void);
};

Kernel_Source::Kernel_Source(void)
{
   numFiles=7;
   program_size=new size_t[numFiles];
   program_buffer=new const char*[numFiles];
   define_InitSigma();
   define_RungeKutta();
   define_Hexciton();
   define_H_antiHermitian();
   define_HPhononLowTemp();
   define_HPhononDoubleExcitonLowTemp();
   define_multiply_dipole_sigmas();
}


void Kernel_Source::define_InitSigma(void)
{
   program_buffer[0]=new char[440];
   program_size[0]=440;
   program_buffer[0]="\
#ifndef INITSIGMA_CL \n\
#define INITSIGMA_CL \n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission \n\
#ifndef TYPE_DEF_MYREAL \n\
#define TYPE_DEF_MYREAL \n\
typedef double myreal; \n\
typedef double2 myreal2; \n\
#endif \n\
__kernel void execute_InitializeSigma(__global myreal2* INsigma, unsigned long int Nmax) \n\
{ \n\
unsigned long int id=get_global_id(0); \n\
  if(id<Nmax) \n\
  { \n\
     INsigma[id].x=0; \n\
     INsigma[id].y=0; \n\
  } \n\
} \n\
#endif\n";
                  }


                     void Kernel_Source::define_RungeKutta(void)
                  {
                     program_buffer[1]=new char[1141];
                     program_size[1]=1141;
                     program_buffer[1]="\
#ifndef RUNGEKUTTA_CL\n\
#define RUNGEKUTTA_CL\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
__kernel void execute_aXsig1Plussig2(__global myreal2* sigOut, __global myreal2* sigIN1, __global myreal2* sigIN2, myreal a, unsigned long int Nmax)\n\
{\n\
   unsigned long int id=get_global_id(0);\n\
   if(id<Nmax)\n\
   {\n\
      sigOut[id].x=a*sigIN1[id].x+sigIN2[id].x;\n\
      sigOut[id].y=a*sigIN1[id].y+sigIN2[id].y;\n\
   }\n\
}\n\
// computation of sigOut+=1/6*sigIN1+1/3*sigIN2+1/3*sigIN3+1/6*sigIN4\n\
__kernel void execute_ReturnNewSigma(__global myreal2* sigOut, __global myreal2* sigIN1, __global myreal2* sigIN2, __global myreal2* sigIN3, __global myreal2* sigIN4, unsigned long int Nmax)\n\
{\n\
   unsigned long int id=get_global_id(0);\n\
   if(id<Nmax)\n\
   {\n\
      myreal2 update;\n\
      update.x=1./6.*sigIN1[id].x+1./3.*sigIN2[id].x+1./3.*sigIN3[id].x+1./6.*sigIN4[id].x;\n\
      update.y=1./6.*sigIN1[id].y+1./3.*sigIN2[id].y+1./3.*sigIN3[id].y+1./6.*sigIN4[id].y;\n\
      sigOut[id].x+=update.x;\n\
      sigOut[id].y+=update.y;\n\
   }\n\
}\n\
#endif\n\
\n";
}


                  void Kernel_Source::define_Hexciton(void)
{
   program_buffer[2]=new char[11337];
   program_size[2]=11337;
   program_buffer[2]="\
#ifndef HEXCITON_CL\n\
#define HEXCITON_CL\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
int returnMatrixEntry(int i, int j, int Matrixdim)\n\
{\n\
   return i*Matrixdim+j;\n\
}\n\
int MemIDSig(int i, int j, int N)\n\
{\n\
   return i*N+j; // so ist es in klasse Matrix\n\
}\n\
int2 integer_quot_rem(unsigned long int x, unsigned long int a)\n\
{\n\
   int2 n;\n\
   n.x=(int)(x / a);\n\
   unsigned long int v = x - n.x * a;\n\
   if ( v < 0 )\n\
      v += a;\n\
   n.y=(int)v;\n\
   return n;\n\
}\n\
__kernel void execute_Hexciton_loopIn(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Hamiltonian)\n\
{\n\
   int id_tuple=get_global_id(0); // index into Ntuple\n\
   if(id_tuple<Ntuples)\n\
   {\n\
      // find start address of sigma-elements in tuple\n\
      int id_sig=id_tuple*Nsites*Nsites;\n\
      int i,j,k;\n\
      myreal hdt=dt/hbar;\n\
      // perform multiplication (Hamiltonian*INsigmaold[id_tuple]-INsigmaold[id_tuple]*Hamiltonian) with three for loops for each thread\n\
      for(i=0; i<Nsites; i++)\n\
      {\n\
         for(j=0; j<Nsites; j++)\n\
         {\n\
            myreal2 tmp;\n\
            tmp.x=0.0;\n\
            tmp.y=0.0;\n\
            for(k=0; k<Nsites; k++)\n\
            {\n\
               tmp.x+=(Hamiltonian[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].x-INsigmaold[id_sig+i*Nsites+k].x*Hamiltonian[k*Nsites+j].x);\n\
               tmp.x-=(Hamiltonian[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].y-INsigmaold[id_sig+i*Nsites+k].y*Hamiltonian[k*Nsites+j].y);\n\
               tmp.y+=(Hamiltonian[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].y-INsigmaold[id_sig+i*Nsites+k].x*Hamiltonian[k*Nsites+j].y);\n\
               tmp.y+=(Hamiltonian[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].x-INsigmaold[id_sig+i*Nsites+k].y*Hamiltonian[k*Nsites+j].x);\n\
            }\n\
            // multiply with -i dt/hbar\n\
            INsigma[id_sig+i*Nsites+j].x += hdt*tmp.y;\n\
            INsigma[id_sig+i*Nsites+j].y -= hdt*tmp.x;\n\
         }\n\
      }\n\
   }\n\
}\n\
__kernel void execute_Hexciton(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Hamiltonian, __local myreal2* sharedHamSig)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block Nblock=Nsites*Nsites\n\
// lade Hamiltonian in shared Memory\n\
      int idelem=get_local_id(0)+get_local_id(1)*get_local_size(0); //goes from 0 .. (Nsites*Nsites-1)\n\
      int idSigShared=idelem+Nblock;                 //goes form (Nblock=Nsites*Nsites) .. (Nblock+Nsites*Nsites-1)\n\
      int idSig=BlockID*Nblock+idelem;  //memory posiion of actual sigmamatrix (BlockID*Nblock) and its elements(idshared1)\n\
///lade Hamiltonian in shared Memory\n\
      sharedHamSig[idelem].x=Hamiltonian[idelem].x;\n\
      sharedHamSig[idelem].y=Hamiltonian[idelem].y;\n\
// lade altes Sigma in shared Memory\n\
      sharedHamSig[idSigShared].x=INsigmaold[idSig].x;\n\
      sharedHamSig[idSigShared].y=INsigmaold[idSig].y;\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
///definiere SpeicherId Matrixelement INsigma[i,j], Speicherposition idSig von oben gehoert zu\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      myreal2 Sigma_act; //Definiere Sigma_act, das vermeidet haeufigen Speicher zugriff auf Elemente von INsigma\n\
      Sigma_act.x=0.; //INsigma->INsigma+Sigma_act am Ende\n\
      Sigma_act.y=0.;\n\
/// berechne H*sigmaold, sigmaold ist im shared Memory\n\
      int idHam;\n\
      int idsigold;\n\
      for(int k=0; k<Nsites; k++)\n\
      {\n\
         idHam=MemIDSig(i, k, Nsites);          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
         idsigold=MemIDSig(k,j, Nsites)+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
         // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
         Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
         Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
         Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
         Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
      }\n\
      /// berechne-sigmaold*H, sigmaold ist im shared Memory\n\
      for(int k=0; k<Nsites; k++)\n\
      {\n\
         idsigold=MemIDSig(i,k, Nsites)+Nblock; //shared memory position corresponding to sigmaold_ik ,i=threadIdx.x\n\
         idHam=MemIDSig(k, j, Nsites);          // shared memory position corresponding to Ham_kj ,j=threadIdx.y\n\
         Sigma_act.x-=sharedHamSig[idsigold].x*sharedHamSig[idHam].x;\n\
         Sigma_act.x+=sharedHamSig[idsigold].y*sharedHamSig[idHam].y;\n\
         Sigma_act.y-=sharedHamSig[idsigold].x*sharedHamSig[idHam].y;\n\
         Sigma_act.y-=sharedHamSig[idsigold].y*sharedHamSig[idHam].x;\n\
      }\n\
      /// multiply by -i/hbar*dt\n\
      myreal Re=Sigma_act.x;\n\
      myreal Im=Sigma_act.y;\n\
      myreal hdt=dt/hbar;\n\
      Sigma_act.x=Im*hdt;\n\
      Sigma_act.y=-Re*hdt;\n\
      /// Final result\n\
      INsigma[idSig].x+=Sigma_act.x;\n\
      INsigma[idSig].y+=Sigma_act.y;\n\
   }\n\
}\n\
__kernel void execute_HexcitonLargeSystem(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Ham, int Nstrides, __local myreal2* sharedHamSig)\n\
{\n\
//  ab hier werden einzelne threads angesprochen\n\
//  threads innerhalb des selben blocks haben shared-memory\n\
//  uebergabe schared memory, hier wird Hamiltonian geladen\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1);\n\
      int jout;\n\
      int iout;\n\
      int idelemOut;\n\
      int kin;\n\
      int idelemHam;\n\
      int idelemSig;\n\
      int idHam;\n\
      int idsigold;\n\
      int idshared;\n\
      int ishare;\n\
      int jshare;\n\
      double hdt=dt/hbar;\n\
      double Re;\n\
      double Im;\n\
      // here need for loop through all elements\n\
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout\n\
      ishare=get_local_id(1);\n\
      jshare=get_local_id(0);\n\
      idshared=jshare+ishare*get_local_size(0);\n\
      for(int jstride=0; jstride<Nstrides; jstride++)\n\
      {\n\
         for(int istride=0; istride<Nstrides; istride++)\n\
         {\n\
            myreal2 Sigma_act;\n\
            Sigma_act.x=0.;\n\
            Sigma_act.y=0.;\n\
            jout=get_local_id(0)+get_local_size(0)*jstride;\n\
            iout=get_local_id(1)+get_local_size(1)*istride;\n\
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemHam=kin+iout*Nsites; //H_(iout,kin)\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared].x=0.;\n\
                  sharedHamSig[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared].x=Ham[idelemHam].x;\n\
                  sharedHamSig[idshared].y=Ham[idelemHam].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemSig=jout+kin*Nsites+BlockID*Nsites*Nsites;\n\
               if(jout+kin*Nsites>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=0.;\n\
                  sharedHamSig[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=INsigmaold[idelemSig].x;\n\
                  sharedHamSig[idshared+Nblock].y=INsigmaold[idelemSig].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            // Multiply -sigmaold*H\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemSig=kin+iout*Nsites+BlockID*Nsites*Nsites;\n\
               if(kin+iout*Nsites>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared].x=0.;\n\
                  sharedHamSig[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared].x=INsigmaold[idelemSig].x;\n\
                  sharedHamSig[idshared].y=INsigmaold[idelemSig].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemHam=jout+kin*Nsites;\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=0.;\n\
                  sharedHamSig[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=-Ham[idelemHam].x;\n\
                  sharedHamSig[idshared+Nblock].y=-Ham[idelemHam].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            Re=Sigma_act.x;\n\
            Im=Sigma_act.y;\n\
            Sigma_act.x=Im*hdt;\n\
            Sigma_act.y=-Re*hdt;\n\
            /// Final result\n\
            if(jout<Nsites && iout<Nsites)\n\
            {\n\
               //if jout,iout\n\
               INsigma[idelemOut].x+=Sigma_act.x;\n\
               INsigma[idelemOut].y+=Sigma_act.y;\n\
            }\n\
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
         }\n\
      }\n\
   }\n\
}\n\
#endif\n";
                  }



                     void Kernel_Source::define_H_antiHermitian(void)
                  {
                     program_buffer[3]=new char[11324];
                     program_size[3]=11324;
                     program_buffer[3]="\
#ifndef H_ANTIHERMITIAN_CL\n\
#define H_ANTIHERMITIAN_CL\n\
\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
\n\
//int returnMatrixEntry(int i, int j, int Matrixdim)\n\
//{\n\
//   return i*Matrixdim+j;\n\
//\n\
//}\n\
\n\
//int MemIDSig(int i, int j, int N)\n\
//{\n\
//   return i*N+j; // so ist es in klasse Matrix\n\
//}\n\
\n\
//int2 integer_quot_rem(unsigned long int x, unsigned long int a)\n\
//{\n\
//   int2 n;\n\
//   n.x=(int)(x / a);\n\
//   unsigned long int v = x - n.x * a;\n\
//   if ( v < 0 )\n\
//      v += a;\n\
//   n.y=(int)v;\n\
//   return n;\n\
//}\n\
\n\
__kernel void execute_antiHermitian_loopIn(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Gamma)\n\
{\n\
   int id_tuple=get_global_id(0); // index into Ntuple\n\
   if(id_tuple<Ntuples)\n\
   {\n\
      // find start address of sigma-elements in touple\n\
      int id_sig=id_tuple*Nsites*Nsites;\n\
      int i,j,k;\n\
      myreal hdt=dt/hbar;\n\
\n\
      // perform multiplication (Gamma*INsigmaold[id_tuple]-INsigmaold[id_tuple]*Gamma) with three for loops for each thread\n\
      for(i=0; i<Nsites; i++)\n\
      {\n\
         for(j=0; j<Nsites; j++)\n\
         {\n\
            myreal2 tmp;\n\
            tmp.x=0.0;\n\
            tmp.y=0.0;\n\
            for(k=0; k<Nsites; k++)\n\
            {\n\
               tmp.x+=(Gamma[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].x+INsigmaold[id_sig+i*Nsites+k].x*Gamma[k*Nsites+j].x);\n\
               tmp.x-=(Gamma[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].y+INsigmaold[id_sig+i*Nsites+k].y*Gamma[k*Nsites+j].y);\n\
               tmp.y+=(Gamma[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].y+INsigmaold[id_sig+i*Nsites+k].x*Gamma[k*Nsites+j].y);\n\
               tmp.y+=(Gamma[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].x+INsigmaold[id_sig+i*Nsites+k].y*Gamma[k*Nsites+j].x);\n\
            }\n\
            // multiply with -1*dt/hbar\n\
            INsigma[id_sig+i*Nsites+j].x -= hdt*tmp.x;\n\
            INsigma[id_sig+i*Nsites+j].y -= hdt*tmp.y;\n\
         }\n\
      }\n\
   }\n\
}\n\
\n\
\n\
\n\
__kernel void execute_antiHermitian(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Gamma, __local myreal2* sharedHamSig)\n\
{\n\
//\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block Nblock=Nsites*Nsites\n\
// lade Gamma in shared Memory\n\
      int idelem=get_local_id(0)+get_local_id(1)*get_local_size(0); //goes from 0 .. (Nsites*Nsites-1)\n\
      int idSigShared=idelem+Nblock;                 //goes form (Nblock=Nsites*Nsites) .. (Nblock+Nsites*Nsites-1)\n\
      int idSig=BlockID*Nblock+idelem;  //memory posiion of actual sigmamatrix (BlockID*Nblock) and its elements(idshared1)\n\
///lade Gamma in shared Memory\n\
      sharedHamSig[idelem].x=Gamma[idelem].x;\n\
      sharedHamSig[idelem].y=Gamma[idelem].y;\n\
// lade altes Sigma in shared Memory\n\
      sharedHamSig[idSigShared].x=INsigmaold[idSig].x;\n\
      sharedHamSig[idSigShared].y=INsigmaold[idSig].y;\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
///definiere SpeicherId Matrixelement INsigma[i,j], Speicherposition idSig von oben gehoert zu\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      myreal2 Sigma_act; //Definiere Sigma_act, das vermeidet haeufigen Speicher zugriff auf Elemente von INsigma\n\
      Sigma_act.x=0.; //INsigma->INsigma+Sigma_act am Ende\n\
      Sigma_act.y=0.;\n\
/// berechne H*sigmaold, sigmaold ist im shared Memory\n\
      int idHam;\n\
      int idsigold;\n\
      for(int k=0; k<Nsites; k++)\n\
      {\n\
         idHam=MemIDSig(i, k, Nsites);          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
         idsigold=MemIDSig(k,j, Nsites)+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
         // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
         Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
         Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
         Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
         Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
//\n\
      }\n\
      /// berechne +sigmaold*H, sigmaold ist im shared Memory\n\
      for(int k=0; k<Nsites; k++)\n\
      {\n\
         idsigold=MemIDSig(i,k, Nsites)+Nblock; //shared memory position corresponding to sigmaold_ik ,i=threadIdx.x\n\
         idHam=MemIDSig(k, j, Nsites);          // shared memory position corresponding to Ham_kj ,j=threadIdx.y\n\
         Sigma_act.x+=sharedHamSig[idsigold].x*sharedHamSig[idHam].x;\n\
         Sigma_act.x-=sharedHamSig[idsigold].y*sharedHamSig[idHam].y;\n\
         Sigma_act.y+=sharedHamSig[idsigold].x*sharedHamSig[idHam].y;\n\
         Sigma_act.y+=sharedHamSig[idsigold].y*sharedHamSig[idHam].x;\n\
      }\n\
      /// multiply by -1/hbar*dt\n\
      myreal Re=Sigma_act.x;\n\
      myreal Im=Sigma_act.y;\n\
      myreal hdt=dt/hbar;\n\
      Sigma_act.x=-Re*hdt;\n\
      Sigma_act.y=-Im*hdt;\n\
      /// Final result\n\
      INsigma[idSig].x+=Sigma_act.x;\n\
      INsigma[idSig].y+=Sigma_act.y;\n\
   }\n\
}\n\
\n\
\n\
\n\
__kernel void execute_antiHermitianLargeSystem(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Ham, int Nstrides, __local myreal2* sharedHamSig)\n\
{\n\
//  ab hier werden einzelne threads angesprochen\n\
//  threads innerhalb des selben blocks haben shared-memory\n\
//  uebergabe schared memory, hier wird Gamma geladen\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1);\n\
      int jout;\n\
      int iout;\n\
      int idelemOut;\n\
      int kin;\n\
      int idelemHam;\n\
      int idelemSig;\n\
      int idHam;\n\
      int idsigold;\n\
      int idshared;\n\
      int ishare;\n\
      int jshare;\n\
      double hdt=dt/hbar;\n\
      double Re;\n\
      double Im;\n\
      // here need for loop through all elements\n\
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout\n\
      ishare=get_local_id(1);\n\
      jshare=get_local_id(0);\n\
      idshared=jshare+ishare*get_local_size(0);\n\
      for(int jstride=0; jstride<Nstrides; jstride++)\n\
      {\n\
         for(int istride=0; istride<Nstrides; istride++)\n\
         {\n\
            myreal2 Sigma_act;\n\
            Sigma_act.x=0.;\n\
            Sigma_act.y=0.;\n\
            jout=get_local_id(0)+get_local_size(0)*jstride;\n\
            iout=get_local_id(1)+get_local_size(1)*istride;\n\
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemHam=kin+iout*Nsites; //H_(iout,kin)\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared].x=0.;\n\
                  sharedHamSig[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared].x=Ham[idelemHam].x;\n\
                  sharedHamSig[idshared].y=Ham[idelemHam].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemSig=jout+kin*Nsites+BlockID*Nsites*Nsites;\n\
               if(jout+kin*Nsites>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=0.;\n\
                  sharedHamSig[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=INsigmaold[idelemSig].x;\n\
                  sharedHamSig[idshared+Nblock].y=INsigmaold[idelemSig].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            // Multiply +sigmaold*H\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemSig=kin+iout*Nsites+BlockID*Nsites*Nsites;\n\
               if(kin+iout*Nsites>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared].x=0.;\n\
                  sharedHamSig[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared].x=INsigmaold[idelemSig].x;\n\
                  sharedHamSig[idshared].y=INsigmaold[idelemSig].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemHam=jout+kin*Nsites;\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=0.;\n\
                  sharedHamSig[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  sharedHamSig[idshared+Nblock].x=Ham[idelemHam].x;\n\
                  sharedHamSig[idshared+Nblock].y=Ham[idelemHam].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;\n\
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;\n\
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            Re=Sigma_act.x;\n\
            Im=Sigma_act.y;\n\
            Sigma_act.x=-Re*hdt;\n\
            Sigma_act.y=-Im*hdt;\n\
            /// Final result\n\
            if(jout<Nsites && iout<Nsites)\n\
            {\n\
               //if jout,iout\n\
               INsigma[idelemOut].x+=Sigma_act.x;\n\
               INsigma[idelemOut].y+=Sigma_act.y;\n\
            }\n\
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
         }\n\
      }\n\
   }\n\
}\n\
#endif\n";
}



                  void Kernel_Source::define_HPhononLowTemp(void)
{
   program_buffer[4]=new char[13794];
   program_size[4]=13794;
   program_buffer[4]="\
#ifndef HPhononLowTemp_CL\n\
#define HPhononLowTemp_CL\n\
\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
\n\
\n\
__kernel void execute_HPhononLowTemp_Part1(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar,__global myreal* d_Rechi0, __global myreal* d_Imchi0, __global myreal* d_S0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      // transfer actuell tuple to shared memory\n\
      // norm of the actual sigma-tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      // get jsite to evaluate Phi_jsite\n\
      int j=d_listNDLsite[get_local_id(0)];\n\
      int NDLjsite=d_listNDL[get_local_id(0)]; //Number of Drude-Lorentz peaks for site jsite\n\
\n\
      ///Berechnung aller (j,k)-Terme\n\
      myreal2 Signew_jk;\n\
      Signew_jk.x=0.;\n\
      Signew_jk.y=0.;\n\
      ///Berechnung aller (j,k)-Terme\n\
      myreal2 Signew_kj;\n\
      Signew_kj.x=0.;\n\
      Signew_kj.y=0.;\n\
      int k=get_local_id(1);\n\
      int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)\n\
      int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)\n\
      //\n\
      /// Phi Anteile gehen mit sigma^(nj+1)\n\
      int SigIDnjplus1;\n\
      int id;\n\
      //  nur falls Nact<Nmax\n\
      if(Nact<Nmax)\n\
      {\n\
         // sum over Drude-Lorentz peaks\n\
         for(int i=0; i<NDLjsite; i++)\n\
         {\n\
            id=d_Memnplus1start[get_local_id(0)]+i;\n\
            SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
            //\n\
            myreal2 Sigold_jk=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk];\n\
            myreal2 Sigold_kj=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj];\n\
            //\n\
            Signew_jk.y+=1./hbar*Sigold_jk.x;\n\
            Signew_jk.x-=1./hbar*Sigold_jk.y;\n\
            //\n\
            Signew_kj.y-=1./hbar*Sigold_kj.x;\n\
            Signew_kj.x+=1./hbar*Sigold_kj.y;\n\
         }\n\
      }\n\
      /// Theta Anteile gehen mit sigma^(nj+1)\n\
      //\n\
      // Schleife ueber verschidene Lorentz-Drude peaks\n\
      int SigIDnjmin1;\n\
      myreal S0;\n\
      myreal Imchi0;\n\
      myreal Rechi0;\n\
      myreal2 LowTempCor;\n\
      for(int i=0; i<NDLjsite; i++)\n\
      {\n\
         id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1\n\
         SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
         //\n\
         S0=d_S0[id];\n\
         Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
         Rechi0=d_Rechi0[id];\n\
         //\n\
         LowTempCor=d_prefactLowTemp[id];\n\
         //\n\
         int idactTup=BlockID*Ncoupled;\n\
         int nact=AllSigmaTuples[idactTup+id];\n\
         if(nact>0) // Teste of nji>0\n\
         {\n\
            myreal2 Sigold_jk=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk];\n\
            myreal2 Sigold_kj=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj];\n\
\n\
            Signew_jk.y+=S0*nact*Sigold_jk.x;\n\
            Signew_jk.x-=S0*nact*Sigold_jk.y;\n\
            //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
            // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
            Signew_jk.x-=Imchi0*nact*Sigold_jk.x;\n\
            Signew_jk.y-=Imchi0*nact*Sigold_jk.y;\n\
            //Multiplikation mit chi0, fuer realteil chi0\n\
            // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
            Signew_jk.x-=Rechi0*nact*Sigold_jk.y;\n\
            Signew_jk.y+=Rechi0*nact*Sigold_jk.x;\n\
\n\
            Signew_kj.y-=S0*nact*Sigold_kj.x;\n\
            Signew_kj.x+=S0*nact*Sigold_kj.y;\n\
            Signew_kj.x-=Imchi0*nact*Sigold_kj.x;\n\
            Signew_kj.y-=Imchi0*nact*Sigold_kj.y;\n\
            Signew_kj.x-=Rechi0*nact*Sigold_kj.y;\n\
            Signew_kj.y+=Rechi0*nact*Sigold_kj.x;\n\
            ///\n\
            /// Low Temperature Correction Theta Part\n\
            ///\n\
            Signew_jk.x+=nact*(LowTempCor.x*Sigold_jk.x - LowTempCor.y*Sigold_jk.y);\n\
            Signew_jk.y+=nact*(LowTempCor.y*Sigold_jk.x + LowTempCor.x*Sigold_jk.y);\n\
            //\n\
            Signew_kj.x+=-nact*(LowTempCor.x*Sigold_kj.x - LowTempCor.y*Sigold_kj.y);\n\
            Signew_kj.y+=-nact*(LowTempCor.y*Sigold_kj.x + LowTempCor.x*Sigold_kj.y);\n\
         }\n\
      }\n\
      ///\n\
      /// Low Temperature Correction, Terme mit sigma^(nj)\n\
      ///\n\
      if(j != k)\n\
      {\n\
         myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];\n\
         myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];\n\
         //\n\
         myreal2 LowTempCor2=d_prefactLowTemp2[get_local_id(0)];\n\
         //\n\
         Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;\n\
         Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;\n\
         //\n\
         Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;\n\
         Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;\n\
      }\n\
\n\
      //\n\
      // /Aktuallisieren von INsigmanew\n\
      Signew_jk.x*=dt;\n\
      Signew_jk.y*=dt;\n\
      Signew_kj.x*=dt;\n\
      Signew_kj.y*=dt;\n\
      //\n\
      int SigIDnew=BlockID;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;\n\
      // FIXME: TK changed __syncthreads to:\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;\n\
      //\n\
   }\n\
}\n\
\n\
__kernel void execute_HPhononLowTemp_Part2(__global myreal2* INsigmanew, __global myreal2* INsigmaold,  int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      myreal2 prefact;\n\
      prefact.x=d_Tupleprefact[BlockID].x;\n\
      prefact.y=d_Tupleprefact[BlockID].y;\n\
\n\
      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      int sig_ij=i*Nsites+j+sigId;\n\
      myreal2 Signew_jk;\n\
      Signew_jk.x=0.;\n\
      Signew_jk.y=0.;\n\
      Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);\n\
      Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);\n\
      //\n\
      INsigmanew[sig_ij].x+=Signew_jk.x;\n\
      INsigmanew[sig_ij].y+=Signew_jk.y;\n\
\n\
   }\n\
}\n\
\n\
__kernel void execute_HPhononLowTemp_Part1LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar,__global myreal* d_Rechi0, __global myreal* d_Imchi0, __global myreal* d_S0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite,int NsitesCoupled)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
   int BlockID=get_group_id(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      // transfer actuell tuple to shared memory\n\
      // norm of the actual sigma-tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      // get jsite to evaluate Phi_jsite\n\
      for(int jj=0; jj<NsitesCoupled; jj++)\n\
      {\n\
\n\
         int j=d_listNDLsite[jj];\n\
         int NDLjsite=d_listNDL[jj]; //Number of Drude-Lorentz peaks for site jsite\n\
         ///Berechnung aller (j,k)-Terme\n\
         myreal2 Signew_jk;\n\
         Signew_jk.x=0.;\n\
         Signew_jk.y=0.;\n\
         ///Berechnung aller (j,k)-Terme\n\
         myreal2 Signew_kj;\n\
         Signew_kj.x=0.;\n\
         Signew_kj.y=0.;\n\
         int k=get_local_id(0);\n\
         int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)\n\
         int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)\n\
         //\n\
         /// Phi Anteile gehen mit sigma^(nj+1)\n\
         int SigIDnjplus1;\n\
         int id;\n\
         //  nur falls Nact<Nmax\n\
         if(Nact<Nmax)\n\
         {\n\
            // sum over Drude-Lorentz peaks\n\
            for(int i=0; i<NDLjsite; i++)\n\
            {\n\
               id=d_Memnplus1start[jj]+i;\n\
               SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
               //\n\
               myreal2 Sigold_jk=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk];\n\
               myreal2 Sigold_kj=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj];\n\
               //\n\
               Signew_jk.y+=1./hbar*Sigold_jk.x;\n\
               Signew_jk.x-=1./hbar*Sigold_jk.y;\n\
               //\n\
               Signew_kj.y-=1./hbar*Sigold_kj.x;\n\
               Signew_kj.x+=1./hbar*Sigold_kj.y;\n\
            }\n\
         }\n\
         /// Theta Anteile gehen mit sigma^(nj+1)\n\
         //\n\
         // Schleife ueber verschidene Lorentz-Drude peaks\n\
         int SigIDnjmin1;\n\
         myreal S0;\n\
         myreal Imchi0;\n\
         myreal Rechi0;\n\
         myreal2 LowTempCor;\n\
         for(int i=0; i<NDLjsite; i++)\n\
         {\n\
            id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1\n\
            SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
            //\n\
            S0=d_S0[id];\n\
            Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
            Rechi0=d_Rechi0[id];\n\
            //\n\
            LowTempCor=d_prefactLowTemp[id];\n\
            //\n\
            int idactTup=BlockID*Ncoupled;\n\
            int nact=AllSigmaTuples[idactTup+id];\n\
            if(nact>0) // Teste of nji>0\n\
            {\n\
               myreal2 Sigold_jk=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk];\n\
               myreal2 Sigold_kj=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj];\n\
\n\
               Signew_jk.y+=S0*nact*Sigold_jk.x;\n\
               Signew_jk.x-=S0*nact*Sigold_jk.y;\n\
               //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
               // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
               Signew_jk.x-=Imchi0*nact*Sigold_jk.x;\n\
               Signew_jk.y-=Imchi0*nact*Sigold_jk.y;\n\
               //Multiplikation mit chi0, fuer realteil chi0\n\
               // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
               Signew_jk.x-=Rechi0*nact*Sigold_jk.y;\n\
               Signew_jk.y+=Rechi0*nact*Sigold_jk.x;\n\
               Signew_kj.y-=S0*nact*Sigold_kj.x;\n\
               Signew_kj.x+=S0*nact*Sigold_kj.y;\n\
               Signew_kj.x-=Imchi0*nact*Sigold_kj.x;\n\
               Signew_kj.y-=Imchi0*nact*Sigold_kj.y;\n\
               Signew_kj.x-=Rechi0*nact*Sigold_kj.y;\n\
               Signew_kj.y+=Rechi0*nact*Sigold_kj.x;\n\
               ///\n\
               /// Low Temperature Correction Theta Part\n\
               ///\n\
               Signew_jk.x+=nact*(LowTempCor.x*Sigold_jk.x - LowTempCor.y*Sigold_jk.y);\n\
               Signew_jk.y+=nact*(LowTempCor.y*Sigold_jk.x + LowTempCor.x*Sigold_jk.y);\n\
               //\n\
               Signew_kj.x+=-nact*(LowTempCor.x*Sigold_kj.x - LowTempCor.y*Sigold_kj.y);\n\
               Signew_kj.y+=-nact*(LowTempCor.y*Sigold_kj.x + LowTempCor.x*Sigold_kj.y);\n\
            }\n\
         }\n\
         ///\n\
         /// Low Temperature Correction, Terme mit sigma^(nj)\n\
         ///\n\
         if(j != k)\n\
         {\n\
            myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];\n\
            myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];\n\
            //\n\
            myreal2 LowTempCor2=d_prefactLowTemp2[jj];\n\
            //\n\
            Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;\n\
            Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;\n\
            //\n\
            Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;\n\
            Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;\n\
         }\n\
\n\
         //\n\
         // /Aktuallisieren von INsigmanew\n\
         Signew_jk.x*=dt;\n\
         Signew_jk.y*=dt;\n\
         Signew_kj.x*=dt;\n\
         Signew_kj.y*=dt;\n\
         //\n\
         int SigIDnew=BlockID;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;\n\
         // FIXME: TK changed __syncthreads to:\n\
         barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;\n\
         //\n\
      }\n\
   }\n\
}\n\
\n\
__kernel void execute_HPhononLowTemp_Part2LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)\n\
{\n\
   int BlockID=get_group_id(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      myreal2 prefact;\n\
      prefact.x=d_Tupleprefact[BlockID].x;\n\
      prefact.y=d_Tupleprefact[BlockID].y;\n\
\n\
      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix\n\
      int j=get_local_id(0);\n\
      for(int ii=0; ii<Nsites; ii++)\n\
      {\n\
         int sig_ij=ii*Nsites+j+sigId;\n\
         myreal2 Signew_jk;\n\
         Signew_jk.x=0.;\n\
         Signew_jk.y=0.;\n\
         Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);\n\
         Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);\n\
         //\n\
         INsigmanew[sig_ij].x+=Signew_jk.x;\n\
         INsigmanew[sig_ij].y+=Signew_jk.y;\n\
      }\n\
   }\n\
}\n\
\n\
#endif\n\
\n";
                  }



                     void Kernel_Source::define_HPhononDoubleExcitonLowTemp(void)
                  {
                     program_buffer[5]=new char[34651];
                     program_size[5]=34651;
                     program_buffer[5]="\
#ifndef HPhononDoubleExcitonLowTemp_CL\n\
#define HPhononDoubleExcitonLowTemp_CL\n\
\n\
\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
\n\
\n\
int ide(int k, int j, int Nsites1) //returns index in Hamiltonian for double excitation |j,k>\n\
{\n\
   int id=Nsites1+(j-k);\n\
   for(int i=1; i<k; i++)\n\
   {\n\
      id+=Nsites1-i;\n\
   }\n\
   return id;\n\
}\n\
\n\
\n\
__kernel void execute_HPhononDoubleLowTemp_Part1(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start,  __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      // transfer actuell tuple to shared memory\n\
//    int idges=threadIdx.x+threadIdx.y*blockDim.x;\n\
// /*   if (idges<Ncoupled)    //Shared geht hier schief\n\
//    {\n\
//     sharedTuple[idges]=AllSigmaTuples[idactTup+idges];\n\
//    }*/\n\
      // norm of the actual sigma-tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      // get jsite to evaluate Phi_jsite\n\
      int j=d_listNDLsite[get_local_id(0)];\n\
      int NDLjsite=d_listNDL[get_local_id(0)]; //Number of Drude-Lorentz peaks for site jsite\n\
      ///Berechnung aller (j,k)-Terme\n\
      myreal2 Signew_jk;\n\
      Signew_jk.x=0.;\n\
      Signew_jk.y=0.;\n\
      ///Berechnung aller (j,k)-Terme\n\
      myreal2 Signew_kj;\n\
      Signew_kj.x=0.;\n\
      Signew_kj.y=0.;\n\
      int k=get_local_id(1);\n\
      int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)\n\
      int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)\n\
      //\n\
      /// Phi Anteile gehen mit sigma^(nj+1)\n\
      int SigIDnjplus1;\n\
      int id;\n\
      //  nur falls Nact<Nmax\n\
      if(Nact<Nmax)\n\
      {\n\
         // sum over Drude-Lorentz peaks\n\
         for(int i=0; i<NDLjsite; i++)\n\
         {\n\
//       id=0;\n\
            id=d_Memnplus1start[get_local_id(0)]+i;\n\
            SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
            Signew_jk.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;\n\
            Signew_jk.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;\n\
            Signew_kj.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;\n\
            Signew_kj.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;\n\
         }\n\
      }\n\
      /// Theta Anteile gehen mit sigma^(nj-1)\n\
      //\n\
      // Schleife ueber verschidene Lorentz-Drude peaks\n\
      int SigIDnjmin1;\n\
      myreal S0;\n\
      myreal Imchi0;\n\
      myreal Rechi0;\n\
      myreal2 LowTempCor;\n\
      for(int i=0; i<NDLjsite; i++)\n\
      {\n\
//            id=0;\n\
         id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1\n\
         SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
         S0=d_S0[id];\n\
         Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
         Rechi0=d_Rechi0[id];\n\
         LowTempCor=d_prefactLowTemp[id];\n\
         //\n\
         int idactTup=BlockID*Ncoupled;\n\
         int nact=AllSigmaTuples[idactTup+id];\n\
\n\
         if(nact>0) // Teste of nji>0\n\
         {\n\
            Signew_jk.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
            Signew_jk.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
            //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
            // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
            Signew_jk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
            Signew_jk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
            //Multiplikation mit chi0, fuer realteil chi0\n\
            // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
            Signew_jk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
            Signew_jk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
\n\
\n\
            Signew_kj.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
            Signew_kj.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
            Signew_kj.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
            Signew_kj.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
            Signew_kj.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
            Signew_kj.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
\n\
            ///\n\
            /// Low Temperature Correction Theta Part\n\
            ///\n\
            Signew_jk.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);\n\
            Signew_jk.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);\n\
            //\n\
            Signew_kj.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);\n\
            Signew_kj.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);\n\
\n\
         }\n\
      }\n\
      //\n\
      ///\n\
      /// Low Temperature Correction, Terme mit sigma^(nj)\n\
      ///\n\
      if(j != k)\n\
      {\n\
         myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];\n\
         myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];\n\
         //\n\
         myreal2 LowTempCor2=d_prefactLowTemp2[get_local_id(0)];\n\
         //\n\
         Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;\n\
         Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;\n\
         //\n\
         Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;\n\
         Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;\n\
      }\n\
      //\n\
      //\n\
      // /Aktuallisieren von INsigmanew\n\
      Signew_jk.x*=dt;\n\
      Signew_jk.y*=dt;\n\
      Signew_kj.x*=dt;\n\
      Signew_kj.y*=dt;\n\
      //\n\
      int SigIDnew=BlockID;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;\n\
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;\n\
      //\n\
   }\n\
}\n\
\n\
\n\
__kernel void execute_HPhononDouble_DoubleExcitationLowTemp(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Nsites1, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
/// Memory id des actuellen tuples, Ncoupled=NsitesCoupled*NLDpeakspersite\n\
   if(BlockID<Ntuples)\n\
   {\n\
      /// Norm aktuelles tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      //\n\
      /// Berechne Speicher ID der aktuellen Sigmamatrix\n\
      int SigIDnew=BlockID;\n\
\n\
      // get jsite to evaluate Phi_jsite\n\
      int j=d_listNDLsite[get_local_id(0)];\n\
      int NDLjsite=d_listNDL[get_local_id(0)]; //Number of Drude-Lorentz peaks for site jsite\n\
\n\
\n\
      ///Berechnung aller (j,k)beta-Terme\n\
      myreal2 Signew_jkbeta;\n\
      Signew_jkbeta.x=0.;\n\
      Signew_jkbeta.y=0.;\n\
      ///Berechnung aller beta(j,k)-Terme\n\
      myreal2 Signew_betajk;\n\
      Signew_betajk.x=0.;\n\
      Signew_betajk.y=0.;\n\
      //\n\
      /// beta index goes over Nsites\n\
      int beta=get_local_id(1);\n\
      //\n\
      /// BEGIN SUMMATION over k for application V_j rho and rho V_j respectively\n\
      for(int n=1; n<=Nsites1; n++) // beachte 0 eintrag entspricht grundzustand, d.h hier summation geht echt bei 1 los\n\
      {\n\
         int k=n;\n\
         Signew_jkbeta.x=0.;\n\
         Signew_jkbeta.y=0.;\n\
         Signew_betajk.x=0.;\n\
         Signew_betajk.y=0.;\n\
         if( j!= k)\n\
         {\n\
            /// BEGIN j != k\n\
            int jk=-100;\n\
            int Matrixpos_jkbeta=-100;\n\
            int Matrixpos_betajk=-100;\n\
            /// define matrix postions of elements jk,beta and beta,jk\n\
            if (k<j)  // do nothing if k==j no double exited sites\n\
            {\n\
               jk=ide(k,j,Nsites1);\n\
               Matrixpos_jkbeta=jk*Nsites+beta;\n\
               Matrixpos_betajk=beta*Nsites+jk;\n\
            }\n\
            if (k>j)\n\
            {\n\
               jk=ide(j,k,Nsites1);\n\
               Matrixpos_jkbeta=jk*Nsites+beta;\n\
               Matrixpos_betajk=beta*Nsites+jk;\n\
            }\n\
       \n\
            /// Phi Anteile gehen mit sigma^(nj+1)\n\
            int SigIDnjplus1;\n\
            int id;\n\
            //  nur falls Nact<Nmax\n\
            if(Nact<Nmax)\n\
            {\n\
               /// BEGIN Nact< Nmax\n\
               // sum over Drude-Lorentz peaks\n\
               for(int i=0; i<NDLjsite; i++)\n\
               {\n\
                  id=d_Memnplus1start[get_local_id(0)]+i;\n\
                  SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
                  //\n\
                  Signew_jkbeta.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                  Signew_jkbeta.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                  //\n\
                  Signew_betajk.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                  Signew_betajk.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].y;\n\
               }\n\
            } /// END Nact< Nmax\n\
            //\n\
            \n\
       \n\
            /// Theta Anteile gehen mit sigma^(nj-1)\n\
            //\n\
            // Schleife ueber verschidene Lorentz-Drude peaks\n\
            int SigIDnjmin1;\n\
            myreal S0;\n\
            myreal Imchi0;\n\
            myreal Rechi0;\n\
            myreal2 LowTempCor;\n\
            for(int i=0; i<NDLjsite; i++)\n\
            {\n\
               /// Begin loop over Drude-Lorentz peaks\n\
               //\n\
               id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1\n\
               SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
               S0=d_S0[id];\n\
               Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
               Rechi0=d_Rechi0[id];\n\
               LowTempCor=d_prefactLowTemp[id];\n\
               int idactTup=BlockID*Ncoupled;\n\
               int nact=AllSigmaTuples[idactTup+id];\n\
               //\n\
               if(nact>0) // Teste of nji>0\n\
               {\n\
                  Signew_jkbeta.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                  Signew_jkbeta.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                  //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
                  // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
                  Signew_jkbeta.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                  Signew_jkbeta.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                  //Multiplikation mit chi0, fuer realteil chi0\n\
                  // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
                  Signew_jkbeta.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                  Signew_jkbeta.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                  //\n\
                  Signew_betajk.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                  Signew_betajk.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                  Signew_betajk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                  Signew_betajk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                  Signew_betajk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                  Signew_betajk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                  ///\n\
                  /// Low Temperature Correction Theta Part\n\
                  ///\n\
                  Signew_jkbeta.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);\n\
                  Signew_jkbeta.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);\n\
                  //\n\
                  Signew_betajk.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);\n\
                  Signew_betajk.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);\n\
               }\n\
            }/// End loop over Drude-Lorentz peaks\n\
            ///\n\
            \n\
    \n\
            /// Low Temperature Correction Teil2\n\
            ///\n\
            ///\n\
            myreal2 prefactLowTemp2 = d_prefactLowTemp2[get_local_id(0)];\n\
            //prefactLowTemp2.x=1;//d_prefactLowTemp2[get_local_id(0)].x;\n\
            //prefactLowTemp2.y=0;//d_prefactLowTemp2[get_local_id(0)].y;\n\
	if(jk != beta)\n\
    	{\n\
            Signew_jkbeta.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
            Signew_jkbeta.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
            Signew_jkbeta.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
            Signew_jkbeta.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
//        //\n\
            Signew_betajk.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;\n\
            Signew_betajk.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;\n\
            Signew_betajk.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;\n\
            Signew_betajk.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;\n\
  	}\n\
  \n\
            //\n\
            //\n\
            // /Aktuallisieren von INsigmanew\n\
            Signew_jkbeta.x*=dt;\n\
            Signew_jkbeta.y*=dt;\n\
            Signew_betajk.x*=dt;\n\
            Signew_betajk.y*=dt;\n\
            //\n\
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x+=Signew_jkbeta.x;\n\
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y+=Signew_jkbeta.y;\n\
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//            __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x+=Signew_betajk.x;\n\
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y+=Signew_betajk.y;\n\
            //\n\
            \n\
     \n\
         }/// End j != k\n\
      } /// Close loop over n\n\
//\n\
   }\n\
}\n\
\n\
__kernel void execute_HPhononDoubleLowTemp_Part2(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Ncoupled, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact, int Ntuples)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
\n\
      myreal2 prefact;\n\
      prefact.x=d_Tupleprefact[BlockID].x;\n\
      prefact.y=d_Tupleprefact[BlockID].y;\n\
\n\
      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      int sig_ij=i*Nsites+j+sigId;\n\
      myreal2 Signew_jk;\n\
      Signew_jk.x=0.;\n\
      Signew_jk.y=0.;\n\
      Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);\n\
      Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);\n\
      myreal2 Signew_ij_update;\n\
      // NOTE perfomance twice as good if I do first load the element that needs to be updated, then update and then store to global memory\n\
      Signew_ij_update.x=INsigmanew[sig_ij].x;\n\
      Signew_ij_update.y=INsigmanew[sig_ij].y;\n\
      Signew_ij_update.x+=Signew_jk.x;\n\
      Signew_ij_update.y+=Signew_jk.y;\n\
      INsigmanew[sig_ij].x=Signew_ij_update.x;\n\
      INsigmanew[sig_ij].y=Signew_ij_update.y;\n\
//      INsigmanew[sig_ij].x+=Signew_jk.x;\n\
//      INsigmanew[sig_ij].y+=Signew_jk.y;\n\
   }\n\
}\n\
\n\
\n\
__kernel void execute_HPhononDoubleLowTemp_Part1LargeSystem(__global myreal2* INsigmanew, __global  myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite, int NsitesCoupled)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
\n\
   int BlockID=get_group_id(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      // norm of the actual sigma-tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      for(int jj=0; jj<NsitesCoupled; jj++)\n\
      {\n\
         // get jsite to evaluate Phi_jsite\n\
         int j=d_listNDLsite[jj];\n\
         int NDLjsite=d_listNDL[jj]; //Number of Drude-Lorentz peaks for site jsite\n\
\n\
\n\
         ///Berechnung aller (j,k)-Terme\n\
         myreal2 Signew_jk;\n\
         Signew_jk.x=0.;\n\
         Signew_jk.y=0.;\n\
         ///Berechnung aller (j,k)-Terme\n\
         myreal2 Signew_kj;\n\
         Signew_kj.x=0.;\n\
         Signew_kj.y=0.;\n\
         int k=get_local_id(0);\n\
         int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)\n\
         int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)\n\
         //\n\
         /// Phi Anteile gehen mit sigma^(nj+1)\n\
         int SigIDnjplus1;\n\
         int id;\n\
         //  nur falls Nact<Nmax\n\
         if(Nact<Nmax)\n\
         {\n\
            // sum over Drude-Lorentz peaks\n\
            for(int i=0; i<NDLjsite; i++)\n\
            {\n\
//       id=0;\n\
               id=d_Memnplus1start[jj]+i;\n\
               SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
               Signew_jk.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;\n\
               Signew_jk.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;\n\
               Signew_kj.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;\n\
               Signew_kj.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;\n\
            }\n\
         }\n\
         /// Theta Anteile gehen mit sigma^(nj-1)\n\
         //\n\
         // Schleife ueber verschidene Lorentz-Drude peaks\n\
         int SigIDnjmin1;\n\
         myreal S0;\n\
         myreal Imchi0;\n\
         myreal Rechi0;\n\
         myreal2 LowTempCor;\n\
         for(int i=0; i<NDLjsite; i++)\n\
         {\n\
//            id=0;\n\
            id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1\n\
            SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
            S0=d_S0[id];\n\
            Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
            Rechi0=d_Rechi0[id];\n\
            LowTempCor=d_prefactLowTemp[id];\n\
            //\n\
            int idactTup=BlockID*Ncoupled;\n\
            int nact=AllSigmaTuples[idactTup+id];\n\
\n\
            if(nact>0) // Teste of nji>0\n\
            {\n\
               Signew_jk.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
               Signew_jk.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
               //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
               // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
               Signew_jk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
               Signew_jk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
               //Multiplikation mit chi0, fuer realteil chi0\n\
               // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
               Signew_jk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;\n\
               Signew_jk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;\n\
\n\
\n\
               Signew_kj.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
               Signew_kj.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
               Signew_kj.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
               Signew_kj.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
               Signew_kj.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;\n\
               Signew_kj.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;\n\
\n\
               ///\n\
               /// Low Temperature Correction Theta Part\n\
               ///\n\
               Signew_jk.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);\n\
               Signew_jk.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);\n\
               //\n\
               Signew_kj.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);\n\
               Signew_kj.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);\n\
\n\
            }\n\
         }\n\
         //\n\
         ///\n\
         /// Low Temperature Correction, Terme mit sigma^(nj)\n\
         ///\n\
         if(j != k)\n\
         {\n\
            myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];\n\
            myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];\n\
            //\n\
            myreal2 LowTempCor2=d_prefactLowTemp2[jj];\n\
            //\n\
            Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;\n\
            Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;\n\
            //\n\
            Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;\n\
            Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;\n\
         }\n\
         //\n\
         //\n\
         // /Aktuallisieren von INsigmanew\n\
         Signew_jk.x*=dt;\n\
         Signew_jk.y*=dt;\n\
         Signew_kj.x*=dt;\n\
         Signew_kj.y*=dt;\n\
         //\n\
         int SigIDnew=BlockID;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;\n\
         barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//         __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;\n\
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;\n\
         //\n\
      }\n\
   }\n\
}\n\
\n\
\n\
__kernel void execute_HPhononDoubleLowTemp_Part2LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Ncoupled, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)\n\
{\n\
   int BlockID=get_group_id(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
\n\
      myreal2 prefact;\n\
      prefact.x=d_Tupleprefact[BlockID].x;\n\
      prefact.y=d_Tupleprefact[BlockID].y;\n\
\n\
      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix\n\
      int j=get_local_id(0);\n\
      for(int ii=0; ii<Nsites; ii++)\n\
      {\n\
         int sig_ij=ii*Nsites+j+sigId;\n\
         myreal2 Signew_jk;\n\
         Signew_jk.x=0.;\n\
         Signew_jk.y=0.;\n\
         Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);\n\
         Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);\n\
         myreal2 Signew_ij_update;\n\
         // NOTE perfomance twice as good if I do first load the element that needs to be updated, then update and then store to global memory\n\
         Signew_ij_update.x=INsigmanew[sig_ij].x;\n\
         Signew_ij_update.y=INsigmanew[sig_ij].y;\n\
         Signew_ij_update.x+=Signew_jk.x;\n\
         Signew_ij_update.y+=Signew_jk.y;\n\
         INsigmanew[sig_ij].x=Signew_ij_update.x;\n\
         INsigmanew[sig_ij].y=Signew_ij_update.y;\n\
//      INsigmanew[sig_ij].x+=Signew_jk.x;\n\
//      INsigmanew[sig_ij].y+=Signew_jk.y;\n\
      }\n\
   }\n\
}\n\
\n\
__kernel void execute_HPhononDouble_DoubleExcitationLowTempLargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Nsites1, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite, int NsitesCoupled)\n\
{\n\
// get BlockID, which tells you which sigma-matrix you are updating\n\
   int BlockID=get_group_id(0);\n\
/// Memory id des actuellen tuples, Ncoupled=NsitesCoupled*NLDpeakspersite\n\
   if(BlockID<Ntuples)\n\
   {\n\
      /// Norm aktuelles tuple\n\
      int Nact=d_TupleNorm[BlockID];\n\
      //\n\
      /// Berechne Speicher ID der aktuellen Sigmamatrix\n\
      int SigIDnew=BlockID;\n\
      for(int jj=0; jj<NsitesCoupled; jj++)\n\
      {\n\
         // get jsite to evaluate Phi_jsite\n\
         int j=d_listNDLsite[jj];\n\
         int NDLjsite=d_listNDL[jj]; //Number of Drude-Lorentz peaks for site jsite\n\
\n\
\n\
         ///Berechnung aller (j,k)beta-Terme\n\
         myreal2 Signew_jkbeta;\n\
         Signew_jkbeta.x=0.;\n\
         Signew_jkbeta.y=0.;\n\
         ///Berechnung aller beta(j,k)-Terme\n\
         myreal2 Signew_betajk;\n\
         Signew_betajk.x=0.;\n\
         Signew_betajk.y=0.;\n\
         //\n\
         /// beta index goes over Nsites\n\
\n\
         int beta=get_local_id(0);\n\
         //\n\
         /// BEGIN SUMMATION over k for application V_j rho and rho V_j respectively\n\
         for(int n=1; n<=Nsites1; n++) // beachte 0 eintrag entspricht grundzustand, d.h hier summation geht echt bei 1 los\n\
         {\n\
            int k=n;\n\
            Signew_jkbeta.x=0.;\n\
            Signew_jkbeta.y=0.;\n\
            Signew_betajk.x=0.;\n\
            Signew_betajk.y=0.;\n\
            if( j!= k)\n\
            {\n\
               /// BEGIN j != k\n\
               int jk=-100;\n\
               int Matrixpos_jkbeta=-100;\n\
               int Matrixpos_betajk=-100;\n\
               /// define matrix postions of elements jk,beta and beta,jk\n\
               if (k<j)  // do nothing if k==j no double exited sites\n\
               {\n\
                  jk=ide(k,j,Nsites1);\n\
                  Matrixpos_jkbeta=jk*Nsites+beta;\n\
                  Matrixpos_betajk=beta*Nsites+jk;\n\
               }\n\
               if (k>j)\n\
               {\n\
                  jk=ide(j,k,Nsites1);\n\
                  Matrixpos_jkbeta=jk*Nsites+beta;\n\
                  Matrixpos_betajk=beta*Nsites+jk;\n\
               }\n\
               /// Phi Anteile gehen mit sigma^(nj+1)\n\
               int SigIDnjplus1;\n\
               int id;\n\
               //  nur falls Nact<Nmax\n\
               if(Nact<Nmax)\n\
               {\n\
                  /// BEGIN Nact< Nmax\n\
                  // sum over Drude-Lorentz peaks\n\
                  for(int i=0; i<NDLjsite; i++)\n\
                  {\n\
                     id=d_Memnplus1start[jj]+i;\n\
                     SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];\n\
                     //\n\
                     Signew_jkbeta.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                     Signew_jkbeta.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                     //\n\
                     Signew_betajk.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                     Signew_betajk.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                  }\n\
               } /// END Nact< Nmax\n\
               //\n\
               /// Theta Anteile gehen mit sigma^(nj-1)\n\
               //\n\
               // Schleife ueber verschidene Lorentz-Drude peaks\n\
               int SigIDnjmin1;\n\
               myreal S0;\n\
               myreal Imchi0;\n\
               myreal Rechi0;\n\
               myreal2 LowTempCor;\n\
               for(int i=0; i<NDLjsite; i++)\n\
               {\n\
                  /// Begin loop over Drude-Lorentz peaks\n\
                  //\n\
                  id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1\n\
                  SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];\n\
                  S0=d_S0[id];\n\
                  Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1\n\
                  Rechi0=d_Rechi0[id];\n\
                  LowTempCor=d_prefactLowTemp[id];\n\
                  int idactTup=BlockID*Ncoupled;\n\
                  int nact=AllSigmaTuples[idactTup+id];\n\
                  //\n\
                  if(nact>0) // Teste of nji>0\n\
                  {\n\
                     Signew_jkbeta.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                     Signew_jkbeta.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                     //Multiplikation mit chi0, fuer imagiaerteil chi0\n\
                     // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)\n\
                     Signew_jkbeta.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                     Signew_jkbeta.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                     //Multiplikation mit chi0, fuer realteil chi0\n\
                     // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)\n\
                     Signew_jkbeta.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
                     Signew_jkbeta.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
                     //\n\
                     Signew_betajk.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                     Signew_betajk.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                     Signew_betajk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                     Signew_betajk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                     Signew_betajk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;\n\
                     Signew_betajk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;\n\
                     ///\n\
                     /// Low Temperature Correction Theta Part\n\
                     ///\n\
                     Signew_jkbeta.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);\n\
                     Signew_jkbeta.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);\n\
                     //\n\
                     Signew_betajk.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);\n\
                     Signew_betajk.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);\n\
                  }\n\
               }/// End loop over Drude-Lorentz peaks\n\
               ///\n\
               /// Low Temperature Correction Teil2\n\
               ///\n\
               ///\n\
               myreal2 prefactLowTemp2;\n\
               prefactLowTemp2.x=d_prefactLowTemp2[jj].x;\n\
               prefactLowTemp2.y=d_prefactLowTemp2[jj].y;\n\
               if(jk != beta)\n\
    		   {\n\
               Signew_jkbeta.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
               Signew_jkbeta.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
               Signew_jkbeta.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;\n\
               Signew_jkbeta.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;\n\
//        //\n\
               Signew_betajk.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;\n\
               Signew_betajk.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;\n\
               Signew_betajk.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;\n\
               Signew_betajk.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;\n\
               }\n\
               //\n\
               //\n\
               // /Aktuallisieren von INsigmanew\n\
               Signew_jkbeta.x*=dt;\n\
               Signew_jkbeta.y*=dt;\n\
               Signew_betajk.x*=dt;\n\
               Signew_betajk.y*=dt;\n\
               //\n\
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x+=Signew_jkbeta.x;\n\
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y+=Signew_jkbeta.y;\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//               __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j\n\
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x+=Signew_betajk.x;\n\
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y+=Signew_betajk.y;\n\
               //\n\
            }/// End j != k\n\
         } /// Close loop over n\n\
//\n\
      }\n\
   }\n\
}\n\
\n\
#endif\n";
}



                  void Kernel_Source::define_multiply_dipole_sigmas(void)
{
   program_buffer[6]=new char[11986];
   program_size[6]=11986;
   program_buffer[6]="\
#ifndef MULTIPLY_DIPOLE_SIGMAS_CL\n\
#define MULTIPLY_DIPOLE_SIGMAS_CL\n\
\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
\n\
\n\
#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission\n\
\n\
#ifndef TYPE_DEF_MYREAL\n\
#define TYPE_DEF_MYREAL\n\
typedef double myreal;\n\
typedef double2 myreal2;\n\
#endif\n\
\n\
\n\
/// this GPU routine multiplies all sigma matrices with the dipole matrices.\n\
//...Left multiplies mu from left and ..Right multiplies from the right.\n\
\n\
__kernel void execute_SigmaMultiplicationLeft(__global myreal2* sigma, __global myreal2* mu, int Ntuples, __local myreal2* SharedMemMuSigma)\n\
{\n\
//	int BlockID=blockIdx.x+blockIdx.y*gridDim.x;\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
// wo bin ich in der Matrix?\n\
// identifiziere sigma und mu Elemente die zu diesem Matrixelement gehoeren\n\
\n\
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}	A_{ij}=double Entry (i*N+j) where N= Matrix dimension\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      int Nsites=get_local_size(1);\n\
      int idelem=i*Nsites+j;		// idelem=id(Element)\n\
\n\
// lade das Element von mu in den shared memory\n\
      SharedMemMuSigma[idelem].x=mu[idelem].x;\n\
      SharedMemMuSigma[idelem].y=mu[idelem].y;\n\
\n\
// lade das Element von sigma in den shared memory\n\
      SharedMemMuSigma[idelem+Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x;\n\
      SharedMemMuSigma[idelem+Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y;\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//	 multipliziere die Elemente aus dem shared memory und specher das Ergebnis lokal.\n\
      sigma[idelem+BlockID*Nsites*Nsites].x=0;\n\
      sigma[idelem+BlockID*Nsites*Nsites].y=0;\n\
      for (int k=0; k<Nsites; k++)\n\
      {\n\
         int id_ik=i*Nsites+k;\n\
         int id_kj=k*Nsites+j;\n\
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}\n\
         sigma[idelem+BlockID*Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].x - SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].y;	\n\
         sigma[idelem+BlockID*Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].y + SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].x;	\n\
      }\n\
   }\n\
}\n\
\n\
__kernel void execute_SigmaMultiplicationRight(__global myreal2* sigma, __global myreal2* mu, int Ntuples, __local myreal2* SharedMemMuSigma)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
// wo bin ich in der Matrix?\n\
// identifiziere sigma und mu Elemente die zu diesem Matrixelement gehoeren\n\
\n\
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}	A_{ij}=double Entry (i*N+j) where N= Matrix dimension\n\
      int i=get_local_id(1);\n\
      int j=get_local_id(0);\n\
      int Nsites=get_local_size(1);\n\
      int idelem=i*Nsites+j;		// idelem=id(Element)\n\
\n\
// lade das Element von mu in den shared memory\n\
      SharedMemMuSigma[idelem+Nsites*Nsites].x=mu[idelem].x;\n\
      SharedMemMuSigma[idelem+Nsites*Nsites].y=mu[idelem].y;\n\
\n\
// lade das Element von sigma in den shared memory\n\
      SharedMemMuSigma[idelem].x=sigma[idelem+BlockID*Nsites*Nsites].x;\n\
      SharedMemMuSigma[idelem].y=sigma[idelem+BlockID*Nsites*Nsites].y;\n\
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
//	 multipliziere die Elemente aus dem shared memory und specher das Ergebnis lokal.\n\
      sigma[idelem+BlockID*Nsites*Nsites].x=0;\n\
      sigma[idelem+BlockID*Nsites*Nsites].y=0;\n\
      for (int k=0; k<Nsites; k++)\n\
      {\n\
         int id_ik=i*Nsites+k;\n\
         int id_kj=k*Nsites+j;\n\
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}\n\
         sigma[idelem+BlockID*Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].x-SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].y;	\n\
         sigma[idelem+BlockID*Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y+SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].x + SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].y;	\n\
      }\n\
   }\n\
}\n\
\n\
\n\
\n\
\n\
__kernel void execute_SigmaMultiplicationLeftLargeSystem(__global myreal2* sigma, __global myreal2* sigma_help, __global myreal2* mu, int Ntuples, int Nsites, int Nstrides, __local myreal2* SharedMemMuSigma)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block\n\
      int jout;\n\
      int iout;\n\
      int idelemOut;\n\
      int kin;\n\
      int idelemHam;\n\
      int idelemSig;\n\
      int idHam;\n\
      int idsigold;\n\
      int idshared;\n\
      int ishare;\n\
      int jshare;\n\
      // here need for loop through all elements\n\
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout\n\
      ishare=get_local_id(1);\n\
      jshare=get_local_id(0);\n\
      idshared=jshare+ishare*get_local_size(0);\n\
      for(int jstride=0; jstride<Nstrides; jstride++)\n\
      {\n\
         for(int istride=0; istride<Nstrides; istride++)\n\
         {\n\
            double2 Sigma_act;\n\
            Sigma_act.x=0.;\n\
            Sigma_act.y=0.;\n\
            jout=get_local_id(0)+get_local_size(0)*jstride;\n\
            iout=get_local_id(1)+get_local_size(1)*istride;\n\
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemHam=kin+iout*Nsites; //H_(iout,kin)\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  SharedMemMuSigma[idshared].x=0.;\n\
                  SharedMemMuSigma[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  SharedMemMuSigma[idshared].x=mu[idelemHam].x;\n\
                  SharedMemMuSigma[idshared].y=mu[idelemHam].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemSig=jout+kin*Nsites+BlockID*Nsites*Nsites;\n\
               if(jout+kin*Nsites>=Nsites*Nsites)\n\
               {\n\
                  SharedMemMuSigma[idshared+Nblock].x=0.;\n\
                  SharedMemMuSigma[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  SharedMemMuSigma[idshared+Nblock].x=sigma[idelemSig].x;\n\
                  SharedMemMuSigma[idshared+Nblock].y=sigma[idelemSig].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished.\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].x;\n\
                  Sigma_act.x-=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].y;\n\
                  Sigma_act.y+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].y;\n\
                  Sigma_act.y+=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            /// Final result\n\
            if(jout<Nsites && iout<Nsites)\n\
            {\n\
               //if jout,iout\n\
               sigma_help[idelemOut].x=Sigma_act.x;\n\
               sigma_help[idelemOut].y=Sigma_act.y;\n\
            }\n\
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
         }\n\
      }\n\
   }\n\
}\n\
\n\
__kernel void execute_SigmaMultiplicationRightLargeSystem(__global myreal2* sigma, __global myreal2* sigma_help, __global myreal2* mu, int Ntuples, int Nsites, int Nstrides, __local myreal2* SharedMemMuSigma)\n\
{\n\
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);\n\
   if(BlockID<Ntuples)\n\
   {\n\
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block\n\
      int jout;\n\
      int iout;\n\
      int idelemOut;\n\
      int kin;\n\
      int idelemHam;\n\
      int idelemSig;\n\
      int idHam;\n\
      int idsigold;\n\
      int idshared;\n\
      int ishare;\n\
      int jshare;\n\
      // here need for loop through all elements\n\
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout\n\
      ishare=get_local_id(1);\n\
      jshare=get_local_id(0);\n\
      idshared=jshare+ishare*get_local_size(0);\n\
      for(int jstride=0; jstride<Nstrides; jstride++)\n\
      {\n\
         for(int istride=0; istride<Nstrides; istride++)\n\
         {\n\
            double2 Sigma_act;\n\
            Sigma_act.x=0.;\n\
            Sigma_act.y=0.;\n\
            jout=get_local_id(0)+get_local_size(0)*jstride;\n\
            iout=get_local_id(1)+get_local_size(1)*istride;\n\
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;\n\
            // Multiply sigma*mu\n\
            for(int kstride=0; kstride<Nstrides; kstride++)\n\
            {\n\
               kin=jshare+kstride*get_local_size(0);\n\
               idelemSig=kin+iout*Nsites+BlockID*Nsites*Nsites;\n\
               if(kin+iout*Nsites>=Nsites*Nsites)\n\
               {\n\
                  SharedMemMuSigma[idshared].x=0.;\n\
                  SharedMemMuSigma[idshared].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  SharedMemMuSigma[idshared].x=sigma[idelemSig].x;\n\
                  SharedMemMuSigma[idshared].y=sigma[idelemSig].y;\n\
               }\n\
               kin=ishare+kstride*get_local_size(0);\n\
               idelemHam=jout+kin*Nsites;\n\
               if(idelemHam>=Nsites*Nsites)\n\
               {\n\
                  SharedMemMuSigma[idshared+Nblock].x=0.;\n\
                  SharedMemMuSigma[idshared+Nblock].y=0.;\n\
               }\n\
               else\n\
               {\n\
                  SharedMemMuSigma[idshared+Nblock].x=mu[idelemHam].x;\n\
                  SharedMemMuSigma[idshared+Nblock].y=mu[idelemHam].y;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
               for(int k=0; k<get_local_size(0); k++)\n\
               {\n\
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x\n\
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y\n\
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x\n\
                  Sigma_act.x+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].x;\n\
                  Sigma_act.x-=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].y;\n\
                  Sigma_act.y+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].y;\n\
                  Sigma_act.y+=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].x;\n\
               }\n\
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished\n\
            }\n\
            /// Final result\n\
            if(jout<Nsites && iout<Nsites)\n\
            {\n\
               //if jout,iout\n\
               sigma_help[idelemOut].x=Sigma_act.x;\n\
               sigma_help[idelemOut].y=Sigma_act.y;\n\
            }\n\
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);\n\
         }\n\
      }\n\
   }\n\
}\n\
\n\
\n\
\n\
__kernel void execute_sigmahelp_to_sigma(__global myreal2* INsigma, __global myreal2* d_sigma_help, unsigned long int Nmax)\n\
{\n\
   unsigned long int id=get_global_id(0);\n\
   if(id<Nmax)\n\
   {\n\
      INsigma[id].x=d_sigma_help[id].x;//\n\
      INsigma[id].y=d_sigma_help[id].y;//\n\
   }\n\
}\n\
\n\
#endif\n";
                  }


#endif


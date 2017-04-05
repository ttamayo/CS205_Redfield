#ifndef HEXCITON_CL
#define HEXCITON_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif

int returnMatrixEntry(int i, int j, int Matrixdim)
{
   return i*Matrixdim+j;

}

int MemIDSig(int i, int j, int N)
{
   return i*N+j; // so ist es in klasse Matrix
}

int2 integer_quot_rem(unsigned long int x, unsigned long int a)
{
   int2 n;
   n.x=(int)(x / a);
   unsigned long int v = x - n.x * a;
   if ( v < 0 )
      v += a;
   n.y=(int)v;
   return n;
}

__kernel void execute_Hexciton_loopIn(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Hamiltonian)
{
   int id_tuple=get_global_id(0); // index into Ntuple
   if(id_tuple<Ntuples)
   {
      // find start address of sigma-elements in tuple
      int id_sig=id_tuple*Nsites*Nsites;
      int i,j,k;
      myreal hdt=dt/hbar;

      // perform multiplication (Hamiltonian*INsigmaold[id_tuple]-INsigmaold[id_tuple]*Hamiltonian) with three for loops for each thread
      for(i=0; i<Nsites; i++)
      {
         for(j=0; j<Nsites; j++)
         {
            myreal2 tmp;
            tmp.x=0.0;
            tmp.y=0.0;
            for(k=0; k<Nsites; k++)
            {
               tmp.x+=(Hamiltonian[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].x-INsigmaold[id_sig+i*Nsites+k].x*Hamiltonian[k*Nsites+j].x);
               tmp.x-=(Hamiltonian[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].y-INsigmaold[id_sig+i*Nsites+k].y*Hamiltonian[k*Nsites+j].y);
               tmp.y+=(Hamiltonian[i*Nsites+k].x*INsigmaold[id_sig+k*Nsites+j].y-INsigmaold[id_sig+i*Nsites+k].x*Hamiltonian[k*Nsites+j].y);
               tmp.y+=(Hamiltonian[i*Nsites+k].y*INsigmaold[id_sig+k*Nsites+j].x-INsigmaold[id_sig+i*Nsites+k].y*Hamiltonian[k*Nsites+j].x);
            }
            // multiply with -i dt/hbar
            INsigma[id_sig+i*Nsites+j].x += hdt*tmp.y;
            INsigma[id_sig+i*Nsites+j].y -= hdt*tmp.x;
         }
      }
   }
}



__kernel void execute_Hexciton(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Hamiltonian, __local myreal2* sharedHamSig)
{
//
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block Nblock=Nsites*Nsites
// lade Hamiltonian in shared Memory
      int idelem=get_local_id(0)+get_local_id(1)*get_local_size(0); //goes from 0 .. (Nsites*Nsites-1)
      int idSigShared=idelem+Nblock;                 //goes form (Nblock=Nsites*Nsites) .. (Nblock+Nsites*Nsites-1)
      int idSig=BlockID*Nblock+idelem;  //memory posiion of actual sigmamatrix (BlockID*Nblock) and its elements(idshared1)
///lade Hamiltonian in shared Memory
      sharedHamSig[idelem].x=Hamiltonian[idelem].x;
      sharedHamSig[idelem].y=Hamiltonian[idelem].y;
// lade altes Sigma in shared Memory
      sharedHamSig[idSigShared].x=INsigmaold[idSig].x;
      sharedHamSig[idSigShared].y=INsigmaold[idSig].y;
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
///definiere SpeicherId Matrixelement INsigma[i,j], Speicherposition idSig von oben gehoert zu
      int i=get_local_id(1);
      int j=get_local_id(0);
      myreal2 Sigma_act; //Definiere Sigma_act, das vermeidet haeufigen Speicher zugriff auf Elemente von INsigma
      Sigma_act.x=0.; //INsigma->INsigma+Sigma_act am Ende
      Sigma_act.y=0.;
/// berechne H*sigmaold, sigmaold ist im shared Memory
      int idHam;
      int idsigold;
      for(int k=0; k<Nsites; k++)
      {
         idHam=MemIDSig(i, k, Nsites);          //shared memory position corresponding to H_ik ,i=threadIdx.x
         idsigold=MemIDSig(k,j, Nsites)+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y
         // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x
         Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;
         Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;
         Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;
         Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;
//
      }
      /// berechne-sigmaold*H, sigmaold ist im shared Memory
      for(int k=0; k<Nsites; k++)
      {
         idsigold=MemIDSig(i,k, Nsites)+Nblock; //shared memory position corresponding to sigmaold_ik ,i=threadIdx.x
         idHam=MemIDSig(k, j, Nsites);          // shared memory position corresponding to Ham_kj ,j=threadIdx.y
         Sigma_act.x-=sharedHamSig[idsigold].x*sharedHamSig[idHam].x;
         Sigma_act.x+=sharedHamSig[idsigold].y*sharedHamSig[idHam].y;
         Sigma_act.y-=sharedHamSig[idsigold].x*sharedHamSig[idHam].y;
         Sigma_act.y-=sharedHamSig[idsigold].y*sharedHamSig[idHam].x;
      }
      /// multiply by -i/hbar*dt
      myreal Re=Sigma_act.x;
      myreal Im=Sigma_act.y;
      myreal hdt=dt/hbar;
      Sigma_act.x=Im*hdt;
      Sigma_act.y=-Re*hdt;
      /// Final result
      INsigma[idSig].x+=Sigma_act.x;
      INsigma[idSig].y+=Sigma_act.y;
   }
}



__kernel void execute_HexcitonLargeSystem(__global myreal2* INsigma, __global myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, __global myreal2* Ham, int Nstrides, __local myreal2* sharedHamSig)
{
//  ab hier werden einzelne threads angesprochen
//  threads innerhalb des selben blocks haben shared-memory
//  uebergabe schared memory, hier wird Hamiltonian geladen
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      int Nblock=get_local_size(0)*get_local_size(1);
      int jout;
      int iout;
      int idelemOut;
      int kin;
      int idelemHam;
      int idelemSig;
      int idHam;
      int idsigold;
      int idshared;
      int ishare;
      int jshare;
      double hdt=dt/hbar;
      double Re;
      double Im;
      // here need for loop through all elements
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout
      ishare=get_local_id(1);
      jshare=get_local_id(0);
      idshared=jshare+ishare*get_local_size(0);
      for(int jstride=0; jstride<Nstrides; jstride++)
      {
         for(int istride=0; istride<Nstrides; istride++)
         {
            myreal2 Sigma_act;
            Sigma_act.x=0.;
            Sigma_act.y=0.;
            jout=get_local_id(0)+get_local_size(0)*jstride;
            iout=get_local_id(1)+get_local_size(1)*istride;
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;
            for(int kstride=0; kstride<Nstrides; kstride++)
            {
               kin=jshare+kstride*get_local_size(0);
               idelemHam=kin+iout*Nsites; //H_(iout,kin)
               if(idelemHam>=Nsites*Nsites)
               {
                  sharedHamSig[idshared].x=0.;
                  sharedHamSig[idshared].y=0.;
               }
               else
               {
                  sharedHamSig[idshared].x=Ham[idelemHam].x;
                  sharedHamSig[idshared].y=Ham[idelemHam].y;
               }
               kin=ishare+kstride*get_local_size(0);
               idelemSig=jout+kin*Nsites+BlockID*Nsites*Nsites;
               if(jout+kin*Nsites>=Nsites*Nsites)
               {
                  sharedHamSig[idshared+Nblock].x=0.;
                  sharedHamSig[idshared+Nblock].y=0.;
               }
               else
               {
                  sharedHamSig[idshared+Nblock].x=INsigmaold[idelemSig].x;
                  sharedHamSig[idshared+Nblock].y=INsigmaold[idelemSig].y;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
               for(int k=0; k<get_local_size(0); k++)
               {
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished
            }
            // Multiply -sigmaold*H
            for(int kstride=0; kstride<Nstrides; kstride++)
            {
               kin=jshare+kstride*get_local_size(0);
               idelemSig=kin+iout*Nsites+BlockID*Nsites*Nsites;
               if(kin+iout*Nsites>=Nsites*Nsites)
               {
                  sharedHamSig[idshared].x=0.;
                  sharedHamSig[idshared].y=0.;
               }
               else
               {
                  sharedHamSig[idshared].x=INsigmaold[idelemSig].x;
                  sharedHamSig[idshared].y=INsigmaold[idelemSig].y;
               }
               kin=ishare+kstride*get_local_size(0);
               idelemHam=jout+kin*Nsites;
               if(idelemHam>=Nsites*Nsites)
               {
                  sharedHamSig[idshared+Nblock].x=0.;
                  sharedHamSig[idshared+Nblock].y=0.;
               }
               else
               {
                  sharedHamSig[idshared+Nblock].x=-Ham[idelemHam].x;
                  sharedHamSig[idshared+Nblock].y=-Ham[idelemHam].y;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
               for(int k=0; k<get_local_size(0); k++)
               {
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x
                  Sigma_act.x+=sharedHamSig[idHam].x*sharedHamSig[idsigold].x;
                  Sigma_act.x-=sharedHamSig[idHam].y*sharedHamSig[idsigold].y;
                  Sigma_act.y+=sharedHamSig[idHam].x*sharedHamSig[idsigold].y;
                  Sigma_act.y+=sharedHamSig[idHam].y*sharedHamSig[idsigold].x;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished
            }
            Re=Sigma_act.x;
            Im=Sigma_act.y;
            Sigma_act.x=Im*hdt;
            Sigma_act.y=-Re*hdt;
            /// Final result
            if(jout<Nsites && iout<Nsites)
            {
               //if jout,iout
               INsigma[idelemOut].x+=Sigma_act.x;
               INsigma[idelemOut].y+=Sigma_act.y;
            }
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
         }
      }
   }
}
#endif

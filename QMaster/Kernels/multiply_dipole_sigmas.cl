#ifndef MULTIPLY_DIPOLE_SIGMAS_CL
#define MULTIPLY_DIPOLE_SIGMAS_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif


#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif


/// this GPU routine multiplies all sigma matrices with the dipole matrices.
//...Left multiplies mu from left and ..Right multiplies from the right.

__kernel void execute_SigmaMultiplicationLeft(__global myreal2* sigma, __global myreal2* mu, int Ntuples, __local myreal2* SharedMemMuSigma)
{
//	int BlockID=blockIdx.x+blockIdx.y*gridDim.x;
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
// wo bin ich in der Matrix?
// identifiziere sigma und mu Elemente die zu diesem Matrixelement gehoeren

// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}	A_{ij}=double Entry (i*N+j) where N= Matrix dimension
      int i=get_local_id(1);
      int j=get_local_id(0);
      int Nsites=get_local_size(1);
      int idelem=i*Nsites+j;		// idelem=id(Element)

// lade das Element von mu in den shared memory
      SharedMemMuSigma[idelem].x=mu[idelem].x;
      SharedMemMuSigma[idelem].y=mu[idelem].y;

// lade das Element von sigma in den shared memory
      SharedMemMuSigma[idelem+Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x;
      SharedMemMuSigma[idelem+Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y;
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//	 multipliziere die Elemente aus dem shared memory und specher das Ergebnis lokal.
      sigma[idelem+BlockID*Nsites*Nsites].x=0;
      sigma[idelem+BlockID*Nsites*Nsites].y=0;
      for (int k=0; k<Nsites; k++)
      {
         int id_ik=i*Nsites+k;
         int id_kj=k*Nsites+j;
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}
         sigma[idelem+BlockID*Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].x - SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].y;
         sigma[idelem+BlockID*Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].y + SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].x;
      }
   }
}

__kernel void execute_SigmaMultiplicationRight(__global myreal2* sigma, __global myreal2* mu, int Ntuples, __local myreal2* SharedMemMuSigma)
{
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
// wo bin ich in der Matrix?
// identifiziere sigma und mu Elemente die zu diesem Matrixelement gehoeren

// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}	A_{ij}=double Entry (i*N+j) where N= Matrix dimension
      int i=get_local_id(1);
      int j=get_local_id(0);
      int Nsites=get_local_size(1);
      int idelem=i*Nsites+j;		// idelem=id(Element)

// lade das Element von mu in den shared memory
      SharedMemMuSigma[idelem+Nsites*Nsites].x=mu[idelem].x;
      SharedMemMuSigma[idelem+Nsites*Nsites].y=mu[idelem].y;

// lade das Element von sigma in den shared memory
      SharedMemMuSigma[idelem].x=sigma[idelem+BlockID*Nsites*Nsites].x;
      SharedMemMuSigma[idelem].y=sigma[idelem+BlockID*Nsites*Nsites].y;
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//	 multipliziere die Elemente aus dem shared memory und specher das Ergebnis lokal.
      sigma[idelem+BlockID*Nsites*Nsites].x=0;
      sigma[idelem+BlockID*Nsites*Nsites].y=0;
      for (int k=0; k<Nsites; k++)
      {
         int id_ik=i*Nsites+k;
         int id_kj=k*Nsites+j;
// multiply A_{ij}=sum_{k} B_{ik}* C_{kj}
         sigma[idelem+BlockID*Nsites*Nsites].x=sigma[idelem+BlockID*Nsites*Nsites].x+SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].x-SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].y;
         sigma[idelem+BlockID*Nsites*Nsites].y=sigma[idelem+BlockID*Nsites*Nsites].y+SharedMemMuSigma[id_ik].y* SharedMemMuSigma[id_kj+Nsites*Nsites].x + SharedMemMuSigma[id_ik].x* SharedMemMuSigma[id_kj+Nsites*Nsites].y;
      }
   }
}




__kernel void execute_SigmaMultiplicationLeftLargeSystem(__global myreal2* sigma, __global myreal2* sigma_help, __global myreal2* mu, int Ntuples, int Nsites, int Nstrides, __local myreal2* SharedMemMuSigma)
{
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block
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
      // here need for loop through all elements
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout
      ishare=get_local_id(1);
      jshare=get_local_id(0);
      idshared=jshare+ishare*get_local_size(0);
      for(int jstride=0; jstride<Nstrides; jstride++)
      {
         for(int istride=0; istride<Nstrides; istride++)
         {
            double2 Sigma_act;
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
                  SharedMemMuSigma[idshared].x=0.;
                  SharedMemMuSigma[idshared].y=0.;
               }
               else
               {
                  SharedMemMuSigma[idshared].x=mu[idelemHam].x;
                  SharedMemMuSigma[idshared].y=mu[idelemHam].y;
               }
               kin=ishare+kstride*get_local_size(0);
               idelemSig=jout+kin*Nsites+BlockID*Nsites*Nsites;
               if(jout+kin*Nsites>=Nsites*Nsites)
               {
                  SharedMemMuSigma[idshared+Nblock].x=0.;
                  SharedMemMuSigma[idshared+Nblock].y=0.;
               }
               else
               {
                  SharedMemMuSigma[idshared+Nblock].x=sigma[idelemSig].x;
                  SharedMemMuSigma[idshared+Nblock].y=sigma[idelemSig].y;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished.
               for(int k=0; k<get_local_size(0); k++)
               {
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x
                  Sigma_act.x+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].x;
                  Sigma_act.x-=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].y;
                  Sigma_act.y+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].y;
                  Sigma_act.y+=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].x;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished
            }
            /// Final result
            if(jout<Nsites && iout<Nsites)
            {
               //if jout,iout
               sigma_help[idelemOut].x=Sigma_act.x;
               sigma_help[idelemOut].y=Sigma_act.y;
            }
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
         }
      }
   }
}

__kernel void execute_SigmaMultiplicationRightLargeSystem(__global myreal2* sigma, __global myreal2* sigma_help, __global myreal2* mu, int Ntuples, int Nsites, int Nstrides, __local myreal2* SharedMemMuSigma)
{
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      int Nblock=get_local_size(0)*get_local_size(1); //number of threads per block
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
      // here need for loop through all elements
      // where to store the output C_iout,jout=Sum_A,ioutk*Bk,jout
      ishare=get_local_id(1);
      jshare=get_local_id(0);
      idshared=jshare+ishare*get_local_size(0);
      for(int jstride=0; jstride<Nstrides; jstride++)
      {
         for(int istride=0; istride<Nstrides; istride++)
         {
            double2 Sigma_act;
            Sigma_act.x=0.;
            Sigma_act.y=0.;
            jout=get_local_id(0)+get_local_size(0)*jstride;
            iout=get_local_id(1)+get_local_size(1)*istride;
            idelemOut=jout+iout*Nsites+BlockID*Nsites*Nsites;
            // Multiply sigma*mu
            for(int kstride=0; kstride<Nstrides; kstride++)
            {
               kin=jshare+kstride*get_local_size(0);
               idelemSig=kin+iout*Nsites+BlockID*Nsites*Nsites;
               if(kin+iout*Nsites>=Nsites*Nsites)
               {
                  SharedMemMuSigma[idshared].x=0.;
                  SharedMemMuSigma[idshared].y=0.;
               }
               else
               {
                  SharedMemMuSigma[idshared].x=sigma[idelemSig].x;
                  SharedMemMuSigma[idshared].y=sigma[idelemSig].y;
               }
               kin=ishare+kstride*get_local_size(0);
               idelemHam=jout+kin*Nsites;
               if(idelemHam>=Nsites*Nsites)
               {
                  SharedMemMuSigma[idshared+Nblock].x=0.;
                  SharedMemMuSigma[idshared+Nblock].y=0.;
               }
               else
               {
                  SharedMemMuSigma[idshared+Nblock].x=mu[idelemHam].x;
                  SharedMemMuSigma[idshared+Nblock].y=mu[idelemHam].y;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
               for(int k=0; k<get_local_size(0); k++)
               {
                  idHam=MemIDSig(ishare, k, get_local_size(0));          //shared memory position corresponding to H_ik ,i=threadIdx.x
                  idsigold=MemIDSig(k,jshare, get_local_size(0))+Nblock; //shared memory position corresponding to sigmaold_kj ,j=threadIdx.y
                  // A.x -> realteil, A.y -> imaginaerteil, bei Multiplikation: C=A*B-> C.x=A.x*B.x-A.y*B.y, C.y=A.x*B.y+A.y*B.x
                  Sigma_act.x+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].x;
                  Sigma_act.x-=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].y;
                  Sigma_act.y+=SharedMemMuSigma[idHam].x*SharedMemMuSigma[idsigold].y;
                  Sigma_act.y+=SharedMemMuSigma[idHam].y*SharedMemMuSigma[idsigold].x;
               }
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);// Wait until first part finished. otherwise later on change shared memory before previous calculations finished
            }
            /// Final result
            if(jout<Nsites && iout<Nsites)
            {
               //if jout,iout
               sigma_help[idelemOut].x=Sigma_act.x;
               sigma_help[idelemOut].y=Sigma_act.y;
            }
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
         }
      }
   }
}



__kernel void execute_sigmahelp_to_sigma(__global myreal2* INsigma, __global myreal2* d_sigma_help, unsigned long int Nmax)
{
   unsigned long int id=get_global_id(0);
   if(id<Nmax)
   {
      INsigma[id].x=d_sigma_help[id].x;//
      INsigma[id].y=d_sigma_help[id].y;//
   }
}

#endif

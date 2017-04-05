#ifndef HPhononLowTemp_CL
#define HPhononLowTemp_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif


__kernel void execute_HPhononLowTemp_Part1(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar,__global myreal* d_Rechi0, __global myreal* d_Imchi0, __global myreal* d_S0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)
{
// get BlockID, which tells you which sigma-matrix you are updating
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      // transfer actuell tuple to shared memory
      // norm of the actual sigma-tuple
      int Nact=d_TupleNorm[BlockID];
      // get jsite to evaluate Phi_jsite
      int j=d_listNDLsite[get_local_id(0)];
      int NDLjsite=d_listNDL[get_local_id(0)]; //Number of Drude-Lorentz peaks for site jsite

      ///Berechnung aller (j,k)-Terme
      myreal2 Signew_jk;
      Signew_jk.x=0.;
      Signew_jk.y=0.;
      ///Berechnung aller (j,k)-Terme
      myreal2 Signew_kj;
      Signew_kj.x=0.;
      Signew_kj.y=0.;
      int k=get_local_id(1);
      int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)
      int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)
      //
      /// Phi Anteile gehen mit sigma^(nj+1)
      int SigIDnjplus1;
      int id;
      //  nur falls Nact<Nmax
      if(Nact<Nmax)
      {
         // sum over Drude-Lorentz peaks
         for(int i=0; i<NDLjsite; i++)
         {
            id=d_Memnplus1start[get_local_id(0)]+i;
            SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
            //
            myreal2 Sigold_jk=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk];
            myreal2 Sigold_kj=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj];
            //
            Signew_jk.y+=1./hbar*Sigold_jk.x;
            Signew_jk.x-=1./hbar*Sigold_jk.y;
            //
            Signew_kj.y-=1./hbar*Sigold_kj.x;
            Signew_kj.x+=1./hbar*Sigold_kj.y;
         }
      }
      /// Theta Anteile gehen mit sigma^(nj+1)
      //
      // Schleife ueber verschidene Lorentz-Drude peaks
      int SigIDnjmin1;
      myreal S0;
      myreal Imchi0;
      myreal Rechi0;
      myreal2 LowTempCor;
      for(int i=0; i<NDLjsite; i++)
      {
         id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1
         SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
         //
         S0=d_S0[id];
         Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
         Rechi0=d_Rechi0[id];
         //
         LowTempCor=d_prefactLowTemp[id];
         //
         int idactTup=BlockID*Ncoupled;
         int nact=AllSigmaTuples[idactTup+id];
         if(nact>0) // Teste of nji>0
         {
            myreal2 Sigold_jk=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk];
            myreal2 Sigold_kj=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj];

            Signew_jk.y+=S0*nact*Sigold_jk.x;
            Signew_jk.x-=S0*nact*Sigold_jk.y;
            //Multiplikation mit chi0, fuer imagiaerteil chi0
            // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
            Signew_jk.x-=Imchi0*nact*Sigold_jk.x;
            Signew_jk.y-=Imchi0*nact*Sigold_jk.y;
            //Multiplikation mit chi0, fuer realteil chi0
            // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
            Signew_jk.x-=Rechi0*nact*Sigold_jk.y;
            Signew_jk.y+=Rechi0*nact*Sigold_jk.x;

            Signew_kj.y-=S0*nact*Sigold_kj.x;
            Signew_kj.x+=S0*nact*Sigold_kj.y;
            Signew_kj.x-=Imchi0*nact*Sigold_kj.x;
            Signew_kj.y-=Imchi0*nact*Sigold_kj.y;
            Signew_kj.x-=Rechi0*nact*Sigold_kj.y;
            Signew_kj.y+=Rechi0*nact*Sigold_kj.x;
            ///
            /// Low Temperature Correction Theta Part
            ///
            Signew_jk.x+=nact*(LowTempCor.x*Sigold_jk.x - LowTempCor.y*Sigold_jk.y);
            Signew_jk.y+=nact*(LowTempCor.y*Sigold_jk.x + LowTempCor.x*Sigold_jk.y);
            //
            Signew_kj.x+=-nact*(LowTempCor.x*Sigold_kj.x - LowTempCor.y*Sigold_kj.y);
            Signew_kj.y+=-nact*(LowTempCor.y*Sigold_kj.x + LowTempCor.x*Sigold_kj.y);
         }
      }
      ///
      /// Low Temperature Correction, Terme mit sigma^(nj)
      ///
      if(j != k)
      {
         myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];
         myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];
         //
         myreal2 LowTempCor2=d_prefactLowTemp2[get_local_id(0)];
         //
         Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;
         Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;
         //
         Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;
         Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;
      }

      //
      // /Aktuallisieren von INsigmanew
      Signew_jk.x*=dt;
      Signew_jk.y*=dt;
      Signew_kj.x*=dt;
      Signew_kj.y*=dt;
      //
      int SigIDnew=BlockID;
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;
      // FIXME: TK changed __syncthreads to:
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;
      //
   }
}

__kernel void execute_HPhononLowTemp_Part2(__global myreal2* INsigmanew, __global myreal2* INsigmaold,  int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)
{
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      myreal2 prefact;
      prefact.x=d_Tupleprefact[BlockID].x;
      prefact.y=d_Tupleprefact[BlockID].y;

      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix
      int i=get_local_id(1);
      int j=get_local_id(0);
      int sig_ij=i*Nsites+j+sigId;
      myreal2 Signew_jk;
      Signew_jk.x=0.;
      Signew_jk.y=0.;
      Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);
      Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);
      //
      INsigmanew[sig_ij].x+=Signew_jk.x;
      INsigmanew[sig_ij].y+=Signew_jk.y;

   }
}

__kernel void execute_HPhononLowTemp_Part1LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar,__global myreal* d_Rechi0, __global myreal* d_Imchi0, __global myreal* d_S0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite,int NsitesCoupled)
{
// get BlockID, which tells you which sigma-matrix you are updating
   int BlockID=get_group_id(0);
   if(BlockID<Ntuples)
   {
      // transfer actuell tuple to shared memory
      // norm of the actual sigma-tuple
      int Nact=d_TupleNorm[BlockID];
      // get jsite to evaluate Phi_jsite
      for(int jj=0; jj<NsitesCoupled; jj++)
      {

         int j=d_listNDLsite[jj];
         int NDLjsite=d_listNDL[jj]; //Number of Drude-Lorentz peaks for site jsite
         ///Berechnung aller (j,k)-Terme
         myreal2 Signew_jk;
         Signew_jk.x=0.;
         Signew_jk.y=0.;
         ///Berechnung aller (j,k)-Terme
         myreal2 Signew_kj;
         Signew_kj.x=0.;
         Signew_kj.y=0.;
         int k=get_local_id(0);
         int Matrixpos_jk=j*Nsites+k; //Speicherelement zu Matrixelement (j,k)
         int Matrixpos_kj=k*Nsites+j; //Speicherelement zu Matrixelement (k,j)
         //
         /// Phi Anteile gehen mit sigma^(nj+1)
         int SigIDnjplus1;
         int id;
         //  nur falls Nact<Nmax
         if(Nact<Nmax)
         {
            // sum over Drude-Lorentz peaks
            for(int i=0; i<NDLjsite; i++)
            {
               id=d_Memnplus1start[jj]+i;
               SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
               //
               myreal2 Sigold_jk=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk];
               myreal2 Sigold_kj=INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj];
               //
               Signew_jk.y+=1./hbar*Sigold_jk.x;
               Signew_jk.x-=1./hbar*Sigold_jk.y;
               //
               Signew_kj.y-=1./hbar*Sigold_kj.x;
               Signew_kj.x+=1./hbar*Sigold_kj.y;
            }
         }
         /// Theta Anteile gehen mit sigma^(nj+1)
         //
         // Schleife ueber verschidene Lorentz-Drude peaks
         int SigIDnjmin1;
         myreal S0;
         myreal Imchi0;
         myreal Rechi0;
         myreal2 LowTempCor;
         for(int i=0; i<NDLjsite; i++)
         {
            id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1
            SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
            //
            S0=d_S0[id];
            Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
            Rechi0=d_Rechi0[id];
            //
            LowTempCor=d_prefactLowTemp[id];
            //
            int idactTup=BlockID*Ncoupled;
            int nact=AllSigmaTuples[idactTup+id];
            if(nact>0) // Teste of nji>0
            {
               myreal2 Sigold_jk=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk];
               myreal2 Sigold_kj=INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj];

               Signew_jk.y+=S0*nact*Sigold_jk.x;
               Signew_jk.x-=S0*nact*Sigold_jk.y;
               //Multiplikation mit chi0, fuer imagiaerteil chi0
               // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
               Signew_jk.x-=Imchi0*nact*Sigold_jk.x;
               Signew_jk.y-=Imchi0*nact*Sigold_jk.y;
               //Multiplikation mit chi0, fuer realteil chi0
               // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
               Signew_jk.x-=Rechi0*nact*Sigold_jk.y;
               Signew_jk.y+=Rechi0*nact*Sigold_jk.x;
               Signew_kj.y-=S0*nact*Sigold_kj.x;
               Signew_kj.x+=S0*nact*Sigold_kj.y;
               Signew_kj.x-=Imchi0*nact*Sigold_kj.x;
               Signew_kj.y-=Imchi0*nact*Sigold_kj.y;
               Signew_kj.x-=Rechi0*nact*Sigold_kj.y;
               Signew_kj.y+=Rechi0*nact*Sigold_kj.x;
               ///
               /// Low Temperature Correction Theta Part
               ///
               Signew_jk.x+=nact*(LowTempCor.x*Sigold_jk.x - LowTempCor.y*Sigold_jk.y);
               Signew_jk.y+=nact*(LowTempCor.y*Sigold_jk.x + LowTempCor.x*Sigold_jk.y);
               //
               Signew_kj.x+=-nact*(LowTempCor.x*Sigold_kj.x - LowTempCor.y*Sigold_kj.y);
               Signew_kj.y+=-nact*(LowTempCor.y*Sigold_kj.x + LowTempCor.x*Sigold_kj.y);
            }
         }
         ///
         /// Low Temperature Correction, Terme mit sigma^(nj)
         ///
         if(j != k)
         {
            myreal2 Sigold_jk=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_jk];
            myreal2 Sigold_kj=INsigmaold[BlockID*Nsites*Nsites+Matrixpos_kj];
            //
            myreal2 LowTempCor2=d_prefactLowTemp2[jj];
            //
            Signew_jk.x+=LowTempCor2.x*Sigold_jk.x - LowTempCor2.y*Sigold_jk.y;
            Signew_jk.y+=LowTempCor2.x*Sigold_jk.y + LowTempCor2.y*Sigold_jk.x;
            //
            Signew_kj.x+=LowTempCor2.x*Sigold_kj.x - LowTempCor2.y*Sigold_kj.y;
            Signew_kj.y+=LowTempCor2.x*Sigold_kj.y + LowTempCor2.y*Sigold_kj.x;
         }

         //
         // /Aktuallisieren von INsigmanew
         Signew_jk.x*=dt;
         Signew_jk.y*=dt;
         Signew_kj.x*=dt;
         Signew_kj.y*=dt;
         //
         int SigIDnew=BlockID;
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].x+=Signew_jk.x;
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jk].y+=Signew_jk.y;
         // FIXME: TK changed __syncthreads to:
         barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;
         //
      }
   }
}

__kernel void execute_HPhononLowTemp_Part2LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)
{
   int BlockID=get_group_id(0);
   if(BlockID<Ntuples)
   {
      myreal2 prefact;
      prefact.x=d_Tupleprefact[BlockID].x;
      prefact.y=d_Tupleprefact[BlockID].y;

      int sigId=BlockID*Nsites*Nsites; //aktuelle Sigmamatrix
      int j=get_local_id(0);
      for(int ii=0; ii<Nsites; ii++)
      {
         int sig_ij=ii*Nsites+j+sigId;
         myreal2 Signew_jk;
         Signew_jk.x=0.;
         Signew_jk.y=0.;
         Signew_jk.x-=dt*(prefact.x*INsigmaold[sig_ij].x-prefact.y*INsigmaold[sig_ij].y);
         Signew_jk.y-=dt*(prefact.x*INsigmaold[sig_ij].y+prefact.y*INsigmaold[sig_ij].x);
         //
         INsigmanew[sig_ij].x+=Signew_jk.x;
         INsigmanew[sig_ij].y+=Signew_jk.y;
      }
   }
}

#endif


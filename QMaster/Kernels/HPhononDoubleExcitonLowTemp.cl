#ifndef HPhononDoubleExcitonLowTemp_CL
#define HPhononDoubleExcitonLowTemp_CL


#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif


int ide(int k, int j, int Nsites1) //returns index in Hamiltonian for double excitation |j,k> if j>k and |k,j> if k>j respectively
{
   int id=Nsites1+(j-k);
   for(int i=1; i<k; i++)
   {
      id+=Nsites1-i;
   }
   return id;
}


__kernel void execute_HPhononDoubleLowTemp_Part1(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start,  __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)
{
// get BlockID, which tells you which sigma-matrix you are updating
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
   if(BlockID<Ntuples)
   {
      // transfer actuell tuple to shared memory
//    int idges=threadIdx.x+threadIdx.y*blockDim.x;
// /*   if (idges<Ncoupled)    //Shared geht hier schief
//    {
//     sharedTuple[idges]=AllSigmaTuples[idactTup+idges];
//    }*/
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
//       id=0;
            id=d_Memnplus1start[get_local_id(0)]+i;
            SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
            Signew_jk.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;
            Signew_jk.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;
            Signew_kj.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;
            Signew_kj.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;
         }
      }
      /// Theta Anteile gehen mit sigma^(nj-1)
      //
      // Schleife ueber verschidene Lorentz-Drude peaks
      int SigIDnjmin1;
      myreal S0;
      myreal Imchi0;
      myreal Rechi0;
      myreal2 LowTempCor;
      for(int i=0; i<NDLjsite; i++)
      {
//            id=0;
         id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1
         SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
         S0=d_S0[id];
         Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
         Rechi0=d_Rechi0[id];
         LowTempCor=d_prefactLowTemp[id];
         //
         int idactTup=BlockID*Ncoupled;
         int nact=AllSigmaTuples[idactTup+id];

         if(nact>0) // Teste of nji>0
         {
            Signew_jk.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;
            Signew_jk.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
            //Multiplikation mit chi0, fuer imagiaerteil chi0
            // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
            Signew_jk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;
            Signew_jk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
            //Multiplikation mit chi0, fuer realteil chi0
            // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
            Signew_jk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
            Signew_jk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;


            Signew_kj.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;
            Signew_kj.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
            Signew_kj.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;
            Signew_kj.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
            Signew_kj.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
            Signew_kj.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;

            ///
            /// Low Temperature Correction Theta Part
            ///
            Signew_jk.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);
            Signew_jk.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);
            //
            Signew_kj.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);
            Signew_kj.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);

         }
      }
      //
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
      barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//      __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;
      INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;
      //
   }
}


__kernel void execute_HPhononDouble_DoubleExcitationLowTemp(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Nsites1, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite)
{
// get BlockID, which tells you which sigma-matrix you are updating
   int BlockID=get_group_id(0)+get_group_id(1)*get_num_groups(0);
/// Memory id des actuellen tuples, Ncoupled=NsitesCoupled*NLDpeakspersite
   if(BlockID<Ntuples)
   {
      /// Norm aktuelles tuple
      int Nact=d_TupleNorm[BlockID];
      //
      /// Berechne Speicher ID der aktuellen Sigmamatrix
      int SigIDnew=BlockID;

      // get jsite to evaluate Phi_jsite
      int j=d_listNDLsite[get_local_id(0)];
      int NDLjsite=d_listNDL[get_local_id(0)]; //Number of Drude-Lorentz peaks for site jsite


      ///Berechnung aller (j,k)beta-Terme
      myreal2 Signew_jkbeta;
      Signew_jkbeta.x=0.;
      Signew_jkbeta.y=0.;
      ///Berechnung aller beta(j,k)-Terme
      myreal2 Signew_betajk;
      Signew_betajk.x=0.;
      Signew_betajk.y=0.;
      //
      /// beta index goes over Nsites
      int beta=get_local_id(1);
      //
      /// BEGIN SUMMATION over k for application V_j rho and rho V_j respectively
      for(int n=1; n<=Nsites1; n++) // beachte 0 eintrag entspricht grundzustand, d.h hier summation geht echt bei 1 los
      {
         int k=n;
         Signew_jkbeta.x=0.;
         Signew_jkbeta.y=0.;
         Signew_betajk.x=0.;
         Signew_betajk.y=0.;
         if( j!= k)
         {
            /// BEGIN j != k
            int jk=-100;
            int Matrixpos_jkbeta=-100;
            int Matrixpos_betajk=-100;
            /// define matrix postions of elements jk,beta and beta,jk
            if (k<j)  // do nothing if k==j no double exited sites
            {
               jk=ide(k,j,Nsites1);
               Matrixpos_jkbeta=jk*Nsites+beta;
               Matrixpos_betajk=beta*Nsites+jk;
            }
            if (k>j)
            {
               jk=ide(j,k,Nsites1);
               Matrixpos_jkbeta=jk*Nsites+beta;
               Matrixpos_betajk=beta*Nsites+jk;
            }
            /// Phi Anteile gehen mit sigma^(nj+1)
            int SigIDnjplus1;
            int id;
            //  nur falls Nact<Nmax
            if(Nact<Nmax)
            {
               /// BEGIN Nact< Nmax
               // sum over Drude-Lorentz peaks
               for(int i=0; i<NDLjsite; i++)
               {
                  id=d_Memnplus1start[get_local_id(0)]+i;
                  SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
                  //
                  Signew_jkbeta.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].x;
                  Signew_jkbeta.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].y;
                  //
                  Signew_betajk.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].x;
                  Signew_betajk.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].y;
               }
            } /// END Nact< Nmax
            //
            /// Theta Anteile gehen mit sigma^(nj-1)
            //
            // Schleife ueber verschidene Lorentz-Drude peaks
            int SigIDnjmin1;
            myreal S0;
            myreal Imchi0;
            myreal Rechi0;
            myreal2 LowTempCor;
            for(int i=0; i<NDLjsite; i++)
            {
               /// Begin loop over Drude-Lorentz peaks
               //
               id=d_Memnplus1start[get_local_id(0)]+i; //Memstart is for nplus1 the same as for nminus1
               SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
               S0=d_S0[id];
               Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
               Rechi0=d_Rechi0[id];
               LowTempCor=d_prefactLowTemp[id];
               int idactTup=BlockID*Ncoupled;
               int nact=AllSigmaTuples[idactTup+id];
               //
               if(nact>0) // Teste of nji>0
               {
                  Signew_jkbeta.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                  Signew_jkbeta.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                  //Multiplikation mit chi0, fuer imagiaerteil chi0
                  // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
                  Signew_jkbeta.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                  Signew_jkbeta.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                  //Multiplikation mit chi0, fuer realteil chi0
                  // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
                  Signew_jkbeta.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                  Signew_jkbeta.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                  //
                  Signew_betajk.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                  Signew_betajk.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                  Signew_betajk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                  Signew_betajk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                  Signew_betajk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                  Signew_betajk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                  ///
                  /// Low Temperature Correction Theta Part
                  ///
                  Signew_jkbeta.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);
                  Signew_jkbeta.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);
                  //
                  Signew_betajk.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);
                  Signew_betajk.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);
               }
            }/// End loop over Drude-Lorentz peaks
            ///
            /// Low Temperature Correction Teil2
            ///
            ///
            myreal2 prefactLowTemp2;
            prefactLowTemp2.x=d_prefactLowTemp2[get_local_id(0)].x;
            prefactLowTemp2.y=d_prefactLowTemp2[get_local_id(0)].y;
            if(jk != beta)
            {
            	Signew_jkbeta.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;
            	Signew_jkbeta.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;
            	Signew_jkbeta.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;
            	Signew_jkbeta.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;
//        //
            	Signew_betajk.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;
            	Signew_betajk.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;
            	Signew_betajk.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;
            	Signew_betajk.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;
			}
            //
            //
            // /Aktuallisieren von INsigmanew
            Signew_jkbeta.x*=dt;
            Signew_jkbeta.y*=dt;
            Signew_betajk.x*=dt;
            Signew_betajk.y*=dt;
            //
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x+=Signew_jkbeta.x;
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y+=Signew_jkbeta.y;
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//            __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x+=Signew_betajk.x;
            INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y+=Signew_betajk.y;
            //
         }/// End j != k
      } /// Close loop over n
//
   }
}

__kernel void execute_HPhononDoubleLowTemp_Part2(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Ncoupled, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact, int Ntuples)
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
      myreal2 Signew_ij_update;
      // NOTE perfomance twice as good if I do first load the element that needs to be updated, then update and then store to global memory
      Signew_ij_update.x=INsigmanew[sig_ij].x;
      Signew_ij_update.y=INsigmanew[sig_ij].y;
      Signew_ij_update.x+=Signew_jk.x;
      Signew_ij_update.y+=Signew_jk.y;
      INsigmanew[sig_ij].x=Signew_ij_update.x;
      INsigmanew[sig_ij].y=Signew_ij_update.y;
//      INsigmanew[sig_ij].x+=Signew_jk.x;
//      INsigmanew[sig_ij].y+=Signew_jk.y;
   }
}


__kernel void execute_HPhononDoubleLowTemp_Part1LargeSystem(__global myreal2* INsigmanew, __global  myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite, int NsitesCoupled)
{
// get BlockID, which tells you which sigma-matrix you are updating

   int BlockID=get_group_id(0);
   if(BlockID<Ntuples)
   {
      // norm of the actual sigma-tuple
      int Nact=d_TupleNorm[BlockID];
      for(int jj=0; jj<NsitesCoupled; jj++)
      {
         // get jsite to evaluate Phi_jsite
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
//       id=0;
               id=d_Memnplus1start[jj]+i;
               SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
               Signew_jk.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].x;
               Signew_jk.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jk].y;
               Signew_kj.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].x;
               Signew_kj.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;//1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_kj].y;
            }
         }
         /// Theta Anteile gehen mit sigma^(nj-1)
         //
         // Schleife ueber verschidene Lorentz-Drude peaks
         int SigIDnjmin1;
         myreal S0;
         myreal Imchi0;
         myreal Rechi0;
         myreal2 LowTempCor;
         for(int i=0; i<NDLjsite; i++)
         {
//            id=0;
            id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1
            SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
            S0=d_S0[id];
            Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
            Rechi0=d_Rechi0[id];
            LowTempCor=d_prefactLowTemp[id];
            //
            int idactTup=BlockID*Ncoupled;
            int nact=AllSigmaTuples[idactTup+id];

            if(nact>0) // Teste of nji>0
            {
               Signew_jk.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;
               Signew_jk.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
               //Multiplikation mit chi0, fuer imagiaerteil chi0
               // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
               Signew_jk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;
               Signew_jk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
               //Multiplikation mit chi0, fuer realteil chi0
               // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
               Signew_jk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y;
               Signew_jk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x;


               Signew_kj.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;
               Signew_kj.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
               Signew_kj.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;
               Signew_kj.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
               Signew_kj.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y;
               Signew_kj.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x;

               ///
               /// Low Temperature Correction Theta Part
               ///
               Signew_jk.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);
               Signew_jk.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jk].y);
               //
               Signew_kj.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);
               Signew_kj.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_kj].y);

            }
         }
         //
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
         barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//         __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].x+=Signew_kj.x;
         INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_kj].y+=Signew_kj.y;
         //
      }
   }
}


__kernel void execute_HPhononDoubleLowTemp_Part2LargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, int Ncoupled, int Nsites,  myreal dt, __global myreal2* d_Tupleprefact,int Ntuples)
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
         myreal2 Signew_ij_update;
         // NOTE perfomance twice as good if I do first load the element that needs to be updated, then update and then store to global memory
         Signew_ij_update.x=INsigmanew[sig_ij].x;
         Signew_ij_update.y=INsigmanew[sig_ij].y;
         Signew_ij_update.x+=Signew_jk.x;
         Signew_ij_update.y+=Signew_jk.y;
         INsigmanew[sig_ij].x=Signew_ij_update.x;
         INsigmanew[sig_ij].y=Signew_ij_update.y;
//      INsigmanew[sig_ij].x+=Signew_jk.x;
//      INsigmanew[sig_ij].y+=Signew_jk.y;
      }
   }
}

__kernel void execute_HPhononDouble_DoubleExcitationLowTempLargeSystem(__global myreal2* INsigmanew, __global myreal2* INsigmaold, __global int* AllSigmaTuples, __global int* d_MemIDnminus1, __global int* d_MemIDnplus1, __global int* d_Memnplus1start, __global int* d_TupleNorm, int Nsites, int Nsites1, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, __global myreal *d_S0, __global myreal *d_Imchi0, __global myreal *d_Rechi0, __global myreal2* d_prefactLowTemp, __global myreal2* d_prefactLowTemp2, __global int* d_listNDL, __global int* d_listNDLsite, int NsitesCoupled)
{
// get BlockID, which tells you which sigma-matrix you are updating
   int BlockID=get_group_id(0);
/// Memory id des actuellen tuples, Ncoupled=NsitesCoupled*NLDpeakspersite
   if(BlockID<Ntuples)
   {
      /// Norm aktuelles tuple
      int Nact=d_TupleNorm[BlockID];
      //
      /// Berechne Speicher ID der aktuellen Sigmamatrix
      int SigIDnew=BlockID;
      for(int jj=0; jj<NsitesCoupled; jj++)
      {
         // get jsite to evaluate Phi_jsite
         int j=d_listNDLsite[jj];
         int NDLjsite=d_listNDL[jj]; //Number of Drude-Lorentz peaks for site jsite


         ///Berechnung aller (j,k)beta-Terme
         myreal2 Signew_jkbeta;
         Signew_jkbeta.x=0.;
         Signew_jkbeta.y=0.;
         ///Berechnung aller beta(j,k)-Terme
         myreal2 Signew_betajk;
         Signew_betajk.x=0.;
         Signew_betajk.y=0.;
         //
         /// beta index goes over Nsites

         int beta=get_local_id(0);
         //
         /// BEGIN SUMMATION over k for application V_j rho and rho V_j respectively
         for(int n=1; n<=Nsites1; n++) // beachte 0 eintrag entspricht grundzustand, d.h hier summation geht echt bei 1 los
         {
            int k=n;
            Signew_jkbeta.x=0.;
            Signew_jkbeta.y=0.;
            Signew_betajk.x=0.;
            Signew_betajk.y=0.;
            if( j!= k)
            {
               /// BEGIN j != k
               int jk=-100;
               int Matrixpos_jkbeta=-100;
               int Matrixpos_betajk=-100;
               /// define matrix postions of elements jk,beta and beta,jk
               if (k<j)  // do nothing if k==j no double exited sites
               {
                  jk=ide(k,j,Nsites1);
                  Matrixpos_jkbeta=jk*Nsites+beta;
                  Matrixpos_betajk=beta*Nsites+jk;
               }
               if (k>j)
               {
                  jk=ide(j,k,Nsites1);
                  Matrixpos_jkbeta=jk*Nsites+beta;
                  Matrixpos_betajk=beta*Nsites+jk;
               }
               /// Phi Anteile gehen mit sigma^(nj+1)
               int SigIDnjplus1;
               int id;
               //  nur falls Nact<Nmax
               if(Nact<Nmax)
               {
                  /// BEGIN Nact< Nmax
                  // sum over Drude-Lorentz peaks
                  for(int i=0; i<NDLjsite; i++)
                  {
                     id=d_Memnplus1start[jj]+i;
                     SigIDnjplus1=d_MemIDnplus1[BlockID*Ncoupled+id];
                     //
                     Signew_jkbeta.y+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].x;
                     Signew_jkbeta.x-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_jkbeta].y;
                     //
                     Signew_betajk.y-=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].x;
                     Signew_betajk.x+=1./hbar*INsigmaold[SigIDnjplus1*Nsites*Nsites+Matrixpos_betajk].y;
                  }
               } /// END Nact< Nmax
               //
               /// Theta Anteile gehen mit sigma^(nj-1)
               //
               // Schleife ueber verschidene Lorentz-Drude peaks
               int SigIDnjmin1;
               myreal S0;
               myreal Imchi0;
               myreal Rechi0;
               myreal2 LowTempCor;
               for(int i=0; i<NDLjsite; i++)
               {
                  /// Begin loop over Drude-Lorentz peaks
                  //
                  id=d_Memnplus1start[jj]+i; //Memstart is for nplus1 the same as for nminus1
                  SigIDnjmin1=d_MemIDnminus1[BlockID*Ncoupled+id];
                  S0=d_S0[id];
                  Imchi0=d_Imchi0[id];  //id geht von 0 bis Ncoupled-1
                  Rechi0=d_Rechi0[id];
                  LowTempCor=d_prefactLowTemp[id];
                  int idactTup=BlockID*Ncoupled;
                  int nact=AllSigmaTuples[idactTup+id];
                  //
                  if(nact>0) // Teste of nji>0
                  {
                     Signew_jkbeta.y+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                     Signew_jkbeta.x-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                     //Multiplikation mit chi0, fuer imagiaerteil chi0
                     // => i*chi0->i*i*Imchi0=-Imchi0 => i*chi0*(a+ib)= -Imchi0*(a+ib)
                     Signew_jkbeta.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                     Signew_jkbeta.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                     //Multiplikation mit chi0, fuer realteil chi0
                     // => i*chi0->i*Rechi0 => i*chi0*(a+ib)= iRechi0*(a+ib)
                     Signew_jkbeta.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y;
                     Signew_jkbeta.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x;
                     //
                     Signew_betajk.y-=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                     Signew_betajk.x+=S0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                     Signew_betajk.x-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                     Signew_betajk.y-=Imchi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                     Signew_betajk.x-=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y;
                     Signew_betajk.y+=Rechi0*nact*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x;
                     ///
                     /// Low Temperature Correction Theta Part
                     ///
                     Signew_jkbeta.x+=nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);
                     Signew_jkbeta.y+=nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_jkbeta].y);
                     //
                     Signew_betajk.x+=-nact*(LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x - LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);
                     Signew_betajk.y+=-nact*(LowTempCor.y*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].x + LowTempCor.x*INsigmaold[SigIDnjmin1*Nsites*Nsites+Matrixpos_betajk].y);
                  }
               }/// End loop over Drude-Lorentz peaks
               ///
               /// Low Temperature Correction Teil2
               ///
               ///
               myreal2 prefactLowTemp2;
               prefactLowTemp2.x=d_prefactLowTemp2[jj].x;
               prefactLowTemp2.y=d_prefactLowTemp2[jj].y;
               if (jk != beta)
               {
               	Signew_jkbeta.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;
               	Signew_jkbeta.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;
               	Signew_jkbeta.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y;
               	Signew_jkbeta.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x;
//        //
              	Signew_betajk.x+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;
               	Signew_betajk.x-=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;
               	Signew_betajk.y+=prefactLowTemp2.x*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y;
               	Signew_betajk.y+=prefactLowTemp2.y*INsigmaold[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x;
				}
               //
               //
               // /Aktuallisieren von INsigmanew
               Signew_jkbeta.x*=dt;
               Signew_jkbeta.y*=dt;
               Signew_betajk.x*=dt;
               Signew_betajk.y*=dt;
               //
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].x+=Signew_jkbeta.x;
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_jkbeta].y+=Signew_jkbeta.y;
               barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
//               __syncthreads(); //wird benoetig um konflikte auszuraeumen fuer k==j
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].x+=Signew_betajk.x;
               INsigmanew[SigIDnew*Nsites*Nsites+Matrixpos_betajk].y+=Signew_betajk.y;
               //
            }/// End j != k
         } /// Close loop over n
//
      }
   }
}

#endif

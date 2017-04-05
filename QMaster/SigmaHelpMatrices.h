#ifndef SIGMAHELPMATRICES_H
#define SIGMAHELPMATRICES_H

#include "Tuple.h"
// // // // // // #include "Matrix.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <map>
// // // // // // #include <complex>
// // // // // //
using namespace std;

class SigmaHelpMatrices
{
public:
   int Nsites; //total sites of the hamiltonian
   int Nmax; //maximal order, hierachically coupled equations are terminated for Sum_i ni >Nmax
   int Ncoupled; //numeber of sites coupling to phonons
   vector<Tuple*> sigmaTuples;
   vector<Tuple*> AllsigmaTuples;  //list of all relevant idices for the sigma matrices, is built after Permutations in sigmaTuples
   int *ListAllsigmaTuples;
   map<string,int> mapSigmaID; //int* sigmaID; //list of memoryID's, sigmaID[id] contains memoryposition of sigma-matrix(n1,n2, ... ,nNsites)
   int *TupleNorm;
   int *MemIDnplus1; //Speicheraddresse der Sigmamatrizen zu n+1
   int *MemIDnminus1; //Speicheraddresse der Sigmamatrizen zu n-1
// // // // // //                       //where id=n1*Nmax^(Nsites-1) + n2*Nmax^(Nsites-2) + ... + (nNsites)
// // // // // //         vector<Matrix*> sigmaMatrices; //list of all used sigma matrices, sigma-matrix(n1,n2, ... ,nNsites), memory acces via sigmaID[id]
public:
   SigmaHelpMatrices(int INNsites, int INNmax, int INNcoupled);
   void translate(Tuple *Intuple, char s[1000]); //is required to create a link between indextuple and sigmamatrix
   void set_sigmaTuples();     //creates index list for sigma matrices, i.e. for 3 sites {{0,0,0},{1,0,0}, {2,0,0}, {1,1,0} ...}
   void set_AllsigmaTuples();    //is created after permutations of sigmaTuples, i.e for 3 sites {{0,0,0},{1,0,0},{0,1,0},{0,0,1}, ... }
   void print_sigmaTuples();  //print vector sigmaTuples;
   void print_AllsigmaTuples();  //print vector AllsigmaTuples;
   long int NumberPermutations(int Ncoupled, int count); //computes Ncoupled!/count! (Ncoupled>=count)
   long int fak(int val);          //computes val!
   long int fak_val1_over_fak_val2(int val1, int val2); //computs val1!/val2!
   void reset_sigmaTuples(); //delete all elements of sigmaTuples, this saves memory
   void initListAllSigma();
// // // // // //         Matrix* returnsigma(Tuple *indices); //Gives sigmaMatrix for certain index-Tupel {n1,n2, ... ,nNsites}
// // // // // //         void assignSigma(Tuple *indices, Matrix *M); //sigmaMatrix for certain index-Tuppel is assignt to Matrix M;
// // // // // //         void assignSigmaEntries(Tuple *indices, Matrix *M); //set all elements of sigmaMatrix for certain index-Tuppel to the elements of Matrix M ;
// // // // // // 	void AddSigmaMatrix(Tuple *indices, Matrix *M); //adds to sigma-matrix for indices Matrix M
// // // // // // 	void ScalarMultiplySigmaMatrix(Tuple *indices, double_complex a); //Multiplyes sigma-matrix for indices with scalar a
// // // // // //
// // // // // // 	double norm(Tuple *indices);
// // // // // // 	void print_sigmaID();
};

SigmaHelpMatrices::SigmaHelpMatrices(int INNsites, int INNmax, int INNcoupled) :
   Nsites(INNsites), Nmax(INNmax), Ncoupled(INNcoupled)
{
   set_sigmaTuples();
   set_AllsigmaTuples();
   reset_sigmaTuples();
   print_AllsigmaTuples();
// // // //   int memsize=pow((Nmax+1),Ncoupled); //+1 because ni=0 is also possible
// // // //   sigmaID=(int*) malloc(sizeof(int) * memsize);
}


void SigmaHelpMatrices::translate(Tuple *Intuple, char s[1000])
{
   sprintf(s,"");
   for(int i=0; i<Intuple->dim; i++)
   {
      sprintf(s,"%s%03d",s,Intuple->entry[i]);
   }
}



void SigmaHelpMatrices::initListAllSigma()
{
   //Memory mapping of sigma matrices
   int count=0; //Counter for memory in sigmaID assocciated with memory postion in sigmaMatrix
/// Create link between sigmamatrises and Index tuple
   if(Nmax>1000)
   {
      fprintf(stdout,"# ERROR SIGMAHELPMATRICES:  Nmax>999\n");
   }
   if(Ncoupled>50000)
   {
      fprintf(stdout,"# ERROR SIGMAHELPMATRICES:  Ncoupled>50000\n");
   }
   char sid[150000];
//   char* sid;
   // ensure that for indextuple ndigits_of_Nmax*Ncoupled fits to character; hier ndigits_of_Nmax restricted to 3
//   sid=new char[3*Ncoupled];
   for(int i=0; i<AllsigmaTuples.size(); i++)
   {
      translate(AllsigmaTuples[i],sid);
//   cout<<sid<<endl;
      mapSigmaID[(string)sid]=count;
//   cout<<mapSigmaID[(string)sid]<<endl;
      count+=1;
// // cout<< "mapSigmaID[" <<id<<"]= "<<mapSigmaID[id]<<endl;
   }
   unsigned long int mem=Ncoupled*AllsigmaTuples.size();
   ListAllsigmaTuples=new int[Ncoupled*AllsigmaTuples.size()];
   count=0;
   for(int i=0; i<AllsigmaTuples.size(); i++)
   {
      for(int j=0; j<Ncoupled; j++)
      {
         ListAllsigmaTuples[count]=AllsigmaTuples[i]->entry[j];
         count+=1;
      }
   }


// Initialisiere zu einzelne Speicher id
// Berechne Norm der Einzelnen Sigmatuples
   TupleNorm=new int[AllsigmaTuples.size()];
   for(int i=0; i<AllsigmaTuples.size(); i++)
   {
      TupleNorm[i]=AllsigmaTuples[i]->sum();
   }

// MemId zu Tuples mit nj+1, Achtung falls Tuplenorm bereits Nmax-> setzte MemId=-1
   MemIDnplus1=new int[mem]; //Speicheraddresse der Sigmamatrizen zu n+1
   for(int i=0; i<AllsigmaTuples.size(); i++)
   {
      for(int j=0; j<Ncoupled; j++) //nj->nj+1
      {
         AllsigmaTuples[i]->entry[j]+=1;
         if(AllsigmaTuples[i]->sum()<=Nmax)
         {
            char sid[1000];
            translate(AllsigmaTuples[i],sid);
            MemIDnplus1[i*Ncoupled+j]=mapSigmaID[(string)sid];
         }
         else MemIDnplus1[i*Ncoupled+j]=-1;
         AllsigmaTuples[i]->entry[j]-=1;
      }
   }

// MemId zu Tuples mit nj-1, Achtung falls nj=0 -> setzte MemId=-1
   MemIDnminus1=new int[mem]; //Speicheraddresse der Sigmamatrizen zu n+1
   for(int i=0; i<AllsigmaTuples.size(); i++)
   {
//     cout<<"# MEMIDnminus1"<<endl;
      for(int j=0; j<Ncoupled; j++) //nj->nj+1
      {
         AllsigmaTuples[i]->entry[j]-=1;
         if(AllsigmaTuples[i]->entry[j]+1>0)
         {
            ///
            char sid[1000];
            translate(AllsigmaTuples[i],sid);
            MemIDnminus1[i*Ncoupled+j] =mapSigmaID[(string)sid];
         }
         else MemIDnminus1[i*Ncoupled+j]=-1;
         AllsigmaTuples[i]->entry[j]+=1;
      }
   }



}


void SigmaHelpMatrices::reset_sigmaTuples()
{
   while(sigmaTuples.size()>1)
   {
      sigmaTuples.pop_back();
   }
}

long int SigmaHelpMatrices::fak(int val)
{
   long int faculty=1;
   if (val==0)
   {
      return 1;
   }
   if (val==1)
   {
      return 1;
   }
   else
   {
      for (int j=1; j<=val; j++)
      {
         faculty*=j;
      }
      return faculty;
   }
}

long int SigmaHelpMatrices::NumberPermutations(int Ncoupled, int count)
{
   long int Npermut=1;
   for(int i=count+1; i<=Ncoupled; i++)
   {
      Npermut*=i;
   }
   return Npermut;
}

long int SigmaHelpMatrices::fak_val1_over_fak_val2(int val1, int val2)
{
   long int faculty=1.;
   if (val2+1>val1)
   {
      return 1;
   }
   else
   {
      for (int j=(val2+1); j<=val1; j++)
      {
         faculty*=j;
      }
      return faculty;
   }
}

void SigmaHelpMatrices::print_sigmaTuples()
{
   for(int i=0; i<sigmaTuples.size(); i++)
   {
      sigmaTuples[i]->print();
   }
}

void SigmaHelpMatrices::print_AllsigmaTuples()
{
//  for(int i=0;i<AllsigmaTuples.size();i++)
//  {
//    AllsigmaTuples[i]->print();
//  }
   fprintf(stdout,"# total number of sigma matrices: %d\n",(int)AllsigmaTuples.size());
   size_t mem=Ncoupled*sizeof(int)*AllsigmaTuples.size();
   fprintf(stdout,"# allocated Device-Memory sigma-tuples: %f MB \n",mem/1024./1024.);
}

void SigmaHelpMatrices::set_sigmaTuples()
{
   //BEGIN set_sigmaTuples
   for(int kk=0; kk<=Nmax; kk++)
   {
      int N=kk;
      const int Ns=Ncoupled;
      vector<Tuple*> HelpTuple;

      int itmin;
      itmin=(N-1)/Ns+1;
      if(N==0)
      {
         itmin=0;
      }
      for(int i=itmin; i<=N; i++)
      {
         Tuple *tuple;
         tuple=new Tuple(Ns);
         tuple->clearEntries();
         tuple->setentry(0,i);
         HelpTuple.push_back(tuple);
      }

      for(int j=0; j<Ns-1; j++)
      {
         int Nmax=HelpTuple.size();

         for(int i=0; i<Nmax; i++)
         {
            for(int k=1; k<Nmax; k++)
            {
               Tuple *tuple;
               tuple=new Tuple(Ns);
               tuple->clearEntries();
               for(int s=0; s<Ns; s++)
               {
                  tuple->entry[s]=HelpTuple[i]->entry[s];
               }
               if(k<=tuple->entry[0+j])
               {
                  tuple->setentry(1+j,k);
                  if(tuple->sum()<=N)
                  {
                     HelpTuple.push_back(tuple);
                  }
               }
            }
         }
      }
      for(int i=0; i<HelpTuple.size(); i++)
      {
         if(HelpTuple[i]->sum()==N)
         {
            sigmaTuples.push_back(HelpTuple[i]);
         }
      }
   }
}//END set_sigmaTuples

void SigmaHelpMatrices::set_AllsigmaTuples()
{
   //BEGIN set_sigmaAllTuples
   for(int ii=0; ii<sigmaTuples.size(); ii++)
   {
      Tuple *tuple;
      tuple=new Tuple(Ncoupled);
      tuple->clearEntries();
      for(int s=0; s<Ncoupled; s++)
      {
         tuple->entry[s]=sigmaTuples[ii]->entry[s];
      }
      //Compute number of possible Permutations
      //Npermut=Nsites! possibilities, introduce factor 1/count! if count numbers are indistinguishable
      int count=1;
      int ilast=0;
      long int Npermut=fak(Ncoupled);
      long int denominator=1;
      for(int i=1; i<Ncoupled; i++)
      {
         if(tuple->entry[i]==tuple->entry[ilast])
         {
            count+=1;
         }
         else
         {
// cout<<"count="<<count<<endl;
            denominator*=fak(count);
//      Npermut=Npermut/fak(count);//NumberPermutations(Ncoupled,count);
//      cout<<ii<<"  hallo1   :"<<count<<endl;
            count=1;
            ilast=i;
         }
      }
      Npermut=fak_val1_over_fak_val2(Ncoupled, count);
      Npermut/=denominator;
//     Npermut=Npermut/fak(count);//fak(Ncoupled-count-1);
// cout<<Npermut<<endl;
      int test[Ncoupled];
      for(int s=0; s<Ncoupled; s++)
      {
         test[s]=tuple->entry[s];
      }

      for(int j=0; j<Npermut; j++)
      {
         next_permutation(test, test + Ncoupled );
         Tuple *tuplehelp;
         tuplehelp=new Tuple(Ncoupled);
         tuplehelp->clearEntries();
         for(int s=0; s<Ncoupled; s++)
         {
            tuplehelp->entry[s]=test[s];
//Generiere hier neues Tupel damit nicht uberschrieben
         }
         AllsigmaTuples.push_back(tuplehelp);
      }
   }
} //END set_sigmaAllTuples


#endif

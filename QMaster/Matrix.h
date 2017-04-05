#ifndef MATRIX_H
#define MATRIX_H

#include "headers.h"
#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>

class Matrix //2-dimensionale Matrix NxN
{
public:
   int N;
   double_complex* entry;
public:
   Matrix(int INN);
   ~Matrix();
   double_complex returnEntry(int i, int j); //gives (i,j)-component of Matrix
   void assignEntry(int i, int j, double_complex a); //assignes (i,j)-component of Matrix with value a
   void assignEntries(Matrix A); //assigne all entries with entries of Matrix A
   void EntryAdd(int i, int j, double_complex a); //add a to  (i,j)-component of Matrix
   void showEntry(int i, int j);
   void clearEntries(); //set all matrix elements to 0
   void print(void);
};

Matrix::Matrix(int INN) :
   N(INN)
{
//Initialise Matrix
   entry=(double_complex*) malloc(sizeof(double_complex)*N*N);
   for(int i=0; i<N*N; i++)
   {
      entry[i]=0.f;
   }
}
Matrix::~Matrix()
{
}


class Vector //2-dimensionale Matrix NxN
{
public:
   int N;
   double_complex* entry;
public:
   Vector(int INN);
   void print(void);
};

Vector::Vector(int INN) :
   N(INN)
{
//Initialise Matrix
   entry=(double_complex*) malloc(sizeof(double_complex)*N);
   for(int i=0; i<N; i++)
   {
      entry[i]=0.;
   }
}

void Vector::print(void)
{
   for(int i=0; i<N; i++)
   {
      printf("(%4e,%4e)", real(entry[i]), imag(entry[i]) );
      printf("\n");
   }
}


void Matrix::print(void)
{
   for(int i=0; i<N; i++)
   {
      for(int j = 0; j < N; j++)
      {
         int id=i*N+j;
//         printf("(%4e,%4e)", real(entry[id]), imag(entry[id]));
         printf("%2e ", real(entry[id]));
      }
      printf("\n");
   }
}

void Matrix::assignEntries(Matrix A)
{
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         int id=i*N+j;
         entry[id]=A.returnEntry(i,j);
      }
   }
}

void Matrix::EntryAdd(int i, int j, double_complex a)
{
   int id=i*N+j;
   entry[id]+=a;
}


void Matrix::clearEntries()
{
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         int id=i*N+j;
         entry[id]=0.;
      }
   }
}

double_complex Matrix::returnEntry(int i, int j)
{
   int id=i*N+j;
   return entry[id];
}

void Matrix::assignEntry(int i, int j, double_complex a)
{
   int id=i*N+j;
   entry[id]=a;
}

void Matrix::showEntry(int i, int j)
{
   int id=i*N+j;
   cout<<"("<<i<<","<<j<<")="<<entry[id]<<endl;
}

void MatrixMuliplication(Matrix *C, Matrix *A, Matrix *B, int N) // C=A*B
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         double_complex c=0.;
         for(int k=0; k<N; k++)
         {
            c+=A->returnEntry(i,k)*B->returnEntry(k,j);
         }
         C->assignEntry(i,j,c);
      }
   }
}

void MatrixAdd(Matrix *C, Matrix *A, Matrix *B, int N)  //C=A+B
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         c=A->returnEntry(i,j)+B->returnEntry(i,j);
         C->assignEntry(i,j,c);
      }
   }
}

void MatrixAdd(Matrix *A, Matrix *B, int N)  //A=A+B
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         c=A->returnEntry(i,j)+B->returnEntry(i,j);
         A->assignEntry(i,j,c);
      }
   }
}

void MatrixSubtract(Matrix *A, Matrix *B, int N)  //A=A-B
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         c=A->returnEntry(i,j)-B->returnEntry(i,j);
         A->assignEntry(i,j,c);
      }
   }
}


void MultiplyScalar(Matrix *A, double_complex a, int N) // A=a*A, a scalar
{
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         double_complex c;
         c=a*A->returnEntry(i,j);
         A->assignEntry(i,j,c);
      }
   }
}

void MultiplyScalar(Matrix *out, Matrix *A, double_complex a, int N) // A=a*A, a scalar
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         c=a*A->returnEntry(i,j);
         out->assignEntry(i,j,c);
      }
   }
}


void Transpose(Matrix *AT, Matrix *A) // AT=Transpose[A]
{
   int N=A->N;
   double_complex c;
   for(int i=0; i<N; i++)
   {
      for(int j=0; j<N; j++)
      {
         c=A->returnEntry(j,i);
         AT->assignEntry(i,j,c);
      }
   }
}


void MatrixVectorMuliplication(Vector *C, Matrix *A, Vector *B, int N) // C=A*B
{
   double_complex c;
   for(int i=0; i<N; i++)
   {
      double_complex c=0.;
      for(int k=0; k<N; k++)
      {
         c+=A->returnEntry(i,k)*B->entry[k];
      }
      C->entry[i]=c;
   }
}


#endif

#ifndef TUPLE_H
#define TUPLE_H

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <complex>
using namespace std;
// Classe Tuple, Definiert ein N-Tupel {n1, n2, n3, ... ,nNmax}
// Funktionen:  sum() bildet Summe ueber alle Eintraege
//    print() gibt Tupel auf Konsole zurueck
class Tuple
{
public:
   int dim;
   int* entry;
public:
   Tuple(int INdim);
//  ~Tuple();
   void print();
   void setentry(int INpos, int val);
   void clearEntries();
   int sum();
};

Tuple::Tuple(int INdim):
   dim(INdim)
{
   entry=(int*) malloc(sizeof(int) * dim);
}

int Tuple::sum()
{
   int s=0;
   for(int i=0; i<dim; i++)
   {
      s+=entry[i];
   }
   return s;
}

void Tuple::clearEntries()
{
   for(int i=0; i <dim; i++)
   {
      entry[i]=0;
   }
}

void Tuple::print()
{
   for(int i = 0; i < dim; i++)
   {
      printf("%4d", entry[i]);
   }
   printf("\n");
}

void Tuple::setentry(int INpos, int val)
{
   entry[INpos]=val;
}


#endif



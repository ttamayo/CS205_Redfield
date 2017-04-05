#ifndef OBSERVABLE_H
#define OBSERVABLE_H


#include <string>
#include <iostream>
using namespace std;

class Observable
{
protected:
   virtual ~Observable();
public:
   virtual void execute(double time, int it, cl_mem d_sigma_System) = 0;
};
Observable::~Observable()
{
}



class ManipulatingObservable
{
protected:
   virtual ~ ManipulatingObservable();
public:
   virtual void execute(double time, int it, cl_mem  d_sigma_System, cl_mem d_sigma_help) = 0;
};
ManipulatingObservable::~ ManipulatingObservable()
{
}



#endif

#ifndef LIOUVILLE_H
#define LIOUVILLE_H

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>
#include "Matrix.h"

using namespace std;

class Liouville
{
protected:
   virtual ~Liouville();
public:
   virtual void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt) = 0; //Berechnet Lroho und gibt output+= Liouv*(INrhoalt,helpMatrices_old) bzw.  helpMatrices_new+=Liouv*(INrhoalt,helpMatrices_old) zurueck
};

Liouville::~Liouville()
{
}



class LiouvilleTime
{
protected:
   virtual ~LiouvilleTime();
public:
   virtual void execute(cl_myreal2* INsigmaold, cl_myreal2* INsigmanew, myreal dt, myreal time) = 0; //Berechnet Lroho und gibt output+= Liouv*(INrhoalt,helpMatrices_old) bzw.  helpMatrices_new+=Liouv*(INrhoalt,helpMatrices_old) zurueck
};

LiouvilleTime::~LiouvilleTime()
{
}

#endif

#ifndef BREAKRULE_H
#define BREAKRULE_H

#include <iostream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <complex>
#include "headers.h"


using namespace std;

class BreakRule
{
protected:
   virtual ~BreakRule();
public:
   virtual void execute(double *breakval, cl_mem d_sigma_System, int it) = 0;
};

BreakRule::~BreakRule()
{
}

#endif

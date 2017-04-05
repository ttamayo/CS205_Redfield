#ifndef INITSIGMA_CL
#define INITSIGMA_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif

__kernel void execute_InitializeSigma(__global myreal2* INsigma, unsigned long int Nmax)
{
   unsigned long int id=get_global_id(0);
   if(id<Nmax)
   {
      INsigma[id].x=0;
      INsigma[id].y=0;
   }
}


#endif

#ifndef RUNGEKUTTA_CL
#define RUNGEKUTTA_CL

#pragma OPENCL EXTENSION cl_khr_fp64 : enable //enable myreal precission

#ifndef TYPE_DEF_MYREAL
#define TYPE_DEF_MYREAL
typedef double myreal;
typedef double2 myreal2;
#endif

__kernel void execute_aXsig1Plussig2(__global myreal2* sigOut, __global myreal2* sigIN1, __global myreal2* sigIN2, myreal a, unsigned long int Nmax)
{
   unsigned long int id=get_global_id(0);
   if(id<Nmax)
   {
      sigOut[id].x=a*sigIN1[id].x+sigIN2[id].x;
      sigOut[id].y=a*sigIN1[id].y+sigIN2[id].y;
   }
}


// computation of sigOut+=1/6*sigIN1+1/3*sigIN2+1/3*sigIN3+1/6*sigIN4
__kernel void execute_ReturnNewSigma(__global myreal2* sigOut, __global myreal2* sigIN1, __global myreal2* sigIN2, __global myreal2* sigIN3, __global myreal2* sigIN4, unsigned long int Nmax)
{
   unsigned long int id=get_global_id(0);
   if(id<Nmax)
   {
      myreal2 update;
      update.x=1./6.*sigIN1[id].x+1./3.*sigIN2[id].x+1./3.*sigIN3[id].x+1./6.*sigIN4[id].x;
      update.y=1./6.*sigIN1[id].y+1./3.*sigIN2[id].y+1./3.*sigIN3[id].y+1./6.*sigIN4[id].y;
      sigOut[id].x+=update.x;
      sigOut[id].y+=update.y;
   }
}



#endif

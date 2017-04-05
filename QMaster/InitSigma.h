#ifndef INITSIGMA_H
#define INITSIGMA_H

#include <iostream>

using namespace std;

void InitializeSigma(cl_myreal2* INsigma, int Ntuples, int Nsites, OpenCL_Init* OpenCLinfo)
{
   int err;
   int workdim=1;
   unsigned long int N=Ntuples*Nsites*Nsites;
   size_t globalWorkSize[]= {N};
   clSetKernelArg(OpenCLinfo->kernel_execute_InitializeSigma,0,sizeof(INsigma), (void*)&INsigma);
   clSetKernelArg(OpenCLinfo->kernel_execute_InitializeSigma,1,sizeof(unsigned long int), (void*)&N);
   cl_event event;
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_InitializeSigma, workdim, 0, globalWorkSize, NULL,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_InitializeSigma failed"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clWaitForEvents(1, &event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_InitializeSigma"<<endl;
      cout<<"Error: clWaitForEvents"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
#ifdef PROFILING
   cl_ulong start_time,end_time;
   size_t return_bytes;
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_START,sizeof(cl_ulong),&start_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_InitializeSigma"<<endl;
      cout<<"Error: clGetEventProfilingInfo (start)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_InitializeSigma"<<endl;
      cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   double run_time =(double)(end_time - start_time);
   fprintf(stderr,"kernel_execute_InitializeSigma: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}




#endif

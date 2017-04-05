#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H
#include <iostream>

void aXsig1Plussig2(cl_myreal2* sigOut, cl_myreal2* sigIN1, cl_myreal2* sigIN2, myreal a, int Ntuples, int Nsites, OpenCL_Init* OpenCLinfo)
{
   unsigned long int N=Ntuples*Nsites*Nsites;
   int err;

   clSetKernelArg(OpenCLinfo->kernel_execute_aXsig1Plussig2,0,sizeof(sigOut), (void*)&sigOut);
   clSetKernelArg(OpenCLinfo->kernel_execute_aXsig1Plussig2,1,sizeof(sigIN1), (void*)&sigIN1);
   clSetKernelArg(OpenCLinfo->kernel_execute_aXsig1Plussig2,2,sizeof(sigIN2), (void*)&sigIN2);
   clSetKernelArg(OpenCLinfo->kernel_execute_aXsig1Plussig2,3,sizeof(myreal), (void*)&a);
   clSetKernelArg(OpenCLinfo->kernel_execute_aXsig1Plussig2,4,sizeof(unsigned long int), (void*)&N);
   cl_event event;
   int workdim=1;
   size_t globalWorkSize[]= {Ntuples*Nsites*Nsites};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_aXsig1Plussig2, workdim, 0, globalWorkSize, NULL,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_aXsig1Plussig2"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clWaitForEvents(1, &event);
   if(err<0)
   {
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
      cout<<"Error: clGetEventProfilingInfo (start)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   double run_time =(double)(end_time - start_time);
   fprintf(stderr,"kernel_execute_aXsig1Plussig2: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}

void ReturnNewSigma(cl_myreal2* sigOut, cl_myreal2* sigIN1, cl_myreal2* sigIN2, cl_myreal2* sigIN3, cl_myreal2* sigIN4, int Ntuples, int Nsites, OpenCL_Init* OpenCLinfo)
{
   unsigned long int N=Ntuples*Nsites*Nsites;
   int err;

   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,0,sizeof(sigOut), (void*)&sigOut);
   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,1,sizeof(sigIN1), (void*)&sigIN1);
   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,2,sizeof(sigIN2), (void*)&sigIN2);
   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,3,sizeof(sigIN3), (void*)&sigIN3);
   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,4,sizeof(sigIN4), (void*)&sigIN4);
   clSetKernelArg(OpenCLinfo->kernel_execute_ReturnNewSigma,5,sizeof(unsigned long int), (void*)&N);
   cl_event event;
   int workdim=1;
   size_t globalWorkSize[]= {Ntuples*Nsites*Nsites};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_ReturnNewSigma, workdim, 0, globalWorkSize,NULL,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_ReturnNewSigma"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clWaitForEvents(1, &event);
   if(err<0)
   {
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
      cout<<"Error: clGetEventProfilingInfo (start)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      exit(1);
   }
   double run_time =(double)(end_time - start_time);
   fprintf(stderr,"kernel_execute_ReturnNewSigma: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


#endif

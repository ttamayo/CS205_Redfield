#ifndef KERNEL_RETURN_WITNESS_H
#define KERNEL_RETURN_WITNESS_H


#include "headers.h"


void Kernel_return_witness(cl_myreal2* d_sigma_System, cl_myreal2* d_sigma_help, cl_myreal2* d_mu, int Nsites, bool left, int Ntuples, OpenCL_Init* OpenCLinfo)
{
   // Multiply with dipoleoperator. Since large system we cannot store the old sigma matrix elements temporary in shared memory. Therefore we use
   // d_sigma_help as container. Later we "copy" the entries of d_sigma_help to d_sigma_system
   // definiere groesse der Thred-Bloecke
   int Nthreadx=16;
   int Nstrides=Nsites/Nthreadx;
   if(Nsites%Nthreadx>0)
   {
      Nstrides+=1;
   }
   int err;
   int workdim=2;
   int dimBlockx=Nthreadx;
   int dimBlocky=Nthreadx;
   size_t localWorkSize[]= {dimBlockx, dimBlocky};
   int ngridx=4*16;
   int ngridy=(Ntuples/ngridx)+1;
   size_t globalWorkSize[]= {ngridx*dimBlockx, ngridy*dimBlocky};
   cl_event event;
   // If multiply from the left
   if (left)
   {
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,0,sizeof(d_sigma_System), (void*)&d_sigma_System);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,1,sizeof(d_sigma_help), (void*)&d_sigma_help);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,2,sizeof(d_mu), (void*)&d_mu);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,3,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,4,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,5,sizeof(int), (void*)&Nstrides);
      // Dynamically allocate local memory (shared memory allocated per workgroup)
      // define shared memory
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem,6,2*Nthreadx*Nthreadx*2*sizeof(myreal),NULL);   // local/shared Memory
      // run kernel
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_SigmaMultiplicationLeftLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: in Kernel_return_witness:"<<endl;
         cout<<"Execute kernel_execute_SigmaMultiplicationLeftLargeSystem failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
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
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
      if(err<0)
      {
         cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      double run_time =(double)(end_time - start_time);
      fprintf(stderr,"kernel_execute_SigmaMultiplicationLeftLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   }
   // If multiply from the right
   if (!left)
   {
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,0,sizeof(d_sigma_System), (void*)&d_sigma_System);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,1,sizeof(d_sigma_help), (void*)&d_sigma_help);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,2,sizeof(d_mu), (void*)&d_mu);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,3,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,4,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,5,sizeof(int), (void*)&Nstrides);
      // Dynamically allocate local memory (shared memory allocated per workgroup)
      // define shared memory
      clSetKernelArg(OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem,6,2*Nthreadx*Nthreadx*2*sizeof(myreal),NULL);   // local/shared Memory
      // run kernel
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_SigmaMultiplicationRightLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: in Kernel_return_witness:"<<endl;
         cout<<"Execute kernel_execute_SigmaMultiplicationRightLargeSystem failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
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
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
      if(err<0)
      {
         cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         cout<<PrintOpenCLerrorMessage(err)<<endl;
         cout.flush();
         exit(1);
      }
      double run_time =(double)(end_time - start_time);
      fprintf(stderr,"kernel_execute_SigmaMultiplicationRightLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   }
   clReleaseEvent(event);
}


#endif

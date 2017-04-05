#ifndef Kernel_anitHermitianRedfield_H
#define Kernel_anitHermitianRedfield_H


// some Notes
//
// BEACHTE KOMPLEXE MULTIPLIKATIONSARITHMETIC RE=RE*RE-IM*IM, IM=RE*IM+IM*RE !!!!!!!
// Hexciton dynamic getestet


void Kernel_antiHermitianRedfield(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Nsites, myreal hbar, myreal dt, cl_myreal2* Gamma, OpenCL_Init* OpenCLinfo)
{
// definiere groesse der Thred-Bloecke
//  int ydim=((Nsites*Nsites)/16)+1*(Nsites*Nsites)%16; //define threadblock with dimension 16 x ydim, try to get muliples of 16, this increases performance of GPU
//  int nblockx=Nsites;//+Nsites%4;
//  int nblocky=Nsites;//+Nsites%2;
   int Ntuples=1;
   int err;
   unsigned long int N=1*Nsites*Nsites;

   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,6,sizeof(Gamma), (void*)&Gamma);
   //
   cl_event event;
   // Dynamically allocate local memory (allocated per workgroup)
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfield,7,2*Nsites*Nsites*2*sizeof(myreal),NULL);   // local/shared Memory
   unsigned long int dimBlockx=Nsites;
   unsigned long int dimBlocky=Nsites;
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=1;
   unsigned long int ngridy=1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitianRedfield, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_antiHermitianRedfield failed"<<endl;
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
   fprintf(stderr,"kernel_execute_antiHermitianRedfield: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


void Kernel_antiHermitianRedfieldLargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Nsites, myreal hbar, myreal dt, cl_myreal2* Gamma, OpenCL_Init* OpenCLinfo)
{
   int err;
   int Ntuples=1; //in Redfield per default Ntuples=1
   //
   cl_event event;
   unsigned long int dimBlockx=16;
   unsigned long int dimBlocky=16;
   if(Nsites>63)
   {
      dimBlockx=32;
      dimBlockx=32;
   }
   int Nstrides=Nsites/dimBlockx;
   if(Nsites%dimBlockx>0)
   {
      Nstrides+=1;
   }
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,6,sizeof(Gamma), (void*)&Gamma);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,7,sizeof(int), (void*)&Nstrides);
   // Dynamically allocate local memory (allocated per workgroup)
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem,8,2*dimBlockx*dimBlocky*2*sizeof(myreal),NULL);   // local/shared Memory
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=1;
   unsigned long int ngridy=1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitianRedfieldLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_antiHermitianRedfieldLargeSystem failed"<<endl;
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
   fprintf(stderr,"kernel_execute_antiHermitianRedfieldLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}



#endif

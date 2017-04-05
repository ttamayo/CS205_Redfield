#ifndef Kernel_anitHermitian_H
#define Kernel_anitHermitian_H


// some Notes
//
// BEACHTE KOMPLEXE MULTIPLIKATIONSARITHMETIC RE=RE*RE-IM*IM, IM=RE*IM+IM*RE !!!!!!!
// Hexciton dynamic getestet


void Kernel_antiHermitian(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, cl_myreal2* Gamma, OpenCL_Init* OpenCLinfo)
{
// definiere groesse der Thred-Bloecke
//  int ydim=((Nsites*Nsites)/16)+1*(Nsites*Nsites)%16; //define threadblock with dimension 16 x ydim, try to get muliples of 16, this increases performance of GPU
//  int nblockx=Nsites;//+Nsites%4;
//  int nblocky=Nsites;//+Nsites%2;

   int err;
   unsigned long int N=Ntuples*Nsites*Nsites;

   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,6,sizeof(Gamma), (void*)&Gamma);
   //
   cl_event event;
   if(strcmp(OpenCLinfo->kernel_type,"LoopIn")==0)
   {
      int workdim=1;
      size_t globalWorkSize[]= {Ntuples};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitian, workdim, 0, globalWorkSize, NULL,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_antiHermitian failed"<<endl;
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
   }

   if(strcmp(OpenCLinfo->kernel_type,"default")==0)
   {
      // Dynamically allocate local memory (allocated per workgroup)
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,7,2*Nsites*Nsites*2*sizeof(myreal),NULL);   // local/shared Memory
      unsigned long int dimBlockx=Nsites;
      unsigned long int dimBlocky=Nsites;
      int workdim=2;
      size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
      unsigned long int ngridx=4*16;
      unsigned long int ngridy=(Ntuples/ngridx)+1;
      size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitian, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_antiHermitian failed"<<endl;
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
   fprintf(stderr,"kernel_execute_antiHermitian: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


void Kernel_antiHermitianLargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, cl_myreal2* Gamma, OpenCL_Init* OpenCLinfo)
{
   int err;
   //
   cl_event event;
   if(strcmp(OpenCLinfo->kernel_type,"LoopIn")==0)
   {
      int workdim=1;
      size_t globalWorkSize[]= {Ntuples};
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,0,sizeof(INsigmanew), (void*)&INsigmanew);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,1,sizeof(INsigmaold), (void*)&INsigmaold);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,2,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,3,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,4,sizeof(myreal), (void*)&hbar);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,5,sizeof(myreal), (void*)&dt);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitian,6,sizeof(Gamma), (void*)&Gamma);
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitian, workdim, 0, globalWorkSize, NULL,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_antiHermitian failed"<<endl;
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
   }

   if(strcmp(OpenCLinfo->kernel_type,"default")==0)
   {
      unsigned long int dimBlockx=16;
      unsigned long int dimBlocky=16;
      int Nstrides=Nsites/dimBlockx;
      if(Nsites%dimBlockx>0)
      {
         Nstrides+=1;
      }
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,0,sizeof(INsigmanew), (void*)&INsigmanew);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,1,sizeof(INsigmaold), (void*)&INsigmaold);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,2,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,3,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,4,sizeof(myreal), (void*)&hbar);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,5,sizeof(myreal), (void*)&dt);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,6,sizeof(Gamma), (void*)&Gamma);
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,7,sizeof(int), (void*)&Nstrides);
      // Dynamically allocate local memory (allocated per workgroup)
      clSetKernelArg(OpenCLinfo->kernel_execute_antiHermitianLargeSystem,8,2*dimBlockx*dimBlocky*2*sizeof(myreal),NULL);   // local/shared Memory
      int workdim=2;
      size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
      unsigned long int ngridx=4*16;
      unsigned long int ngridy=(Ntuples/ngridx)+1;
      size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_antiHermitianLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_antiHermitianLargeSystem failed"<<endl;
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
   fprintf(stderr,"kernel_execute_antiHermitianLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}



#endif

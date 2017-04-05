#ifndef Kernel_HEXCITON_H
#define Kernel_HEXCITON_H


// some Notes
//
// BEACHTE KOMPLEXE MULTIPLIKATIONSARITHMETIC RE=RE*RE-IM*IM, IM=RE*IM+IM*RE !!!!!!!
// Hexciton dynamic getestet


void Kernel_Hexciton(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, cl_myreal2* Ham, OpenCL_Init* OpenCLinfo)
{
// definiere groesse der Thred-Bloecke
//  int ydim=((Nsites*Nsites)/16)+1*(Nsites*Nsites)%16; //define threadblock with dimension 16 x ydim, try to get muliples of 16, this increases performance of GPU
//  int nblockx=Nsites;//+Nsites%4;
//  int nblocky=Nsites;//+Nsites%2;

   int err;
   unsigned long int N=Ntuples*Nsites*Nsites;

   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,6,sizeof(Ham), (void*)&Ham);
   //
   cl_event event;
   if(strcmp(OpenCLinfo->kernel_type,"LoopIn")==0)
   {
      int workdim=1;
      size_t globalWorkSize[]= {Ntuples};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_Hexciton, workdim, 0, globalWorkSize, NULL,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_Hexciton failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
         exit(1);
      }
   }

   if(strcmp(OpenCLinfo->kernel_type,"default")==0)
   {
      // Dynamically allocate local memory (allocated per workgroup)
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,7,2*Nsites*Nsites*2*sizeof(myreal),NULL);   // local/shared Memory
      unsigned long int dimBlockx=Nsites;
      unsigned long int dimBlocky=Nsites;
      int workdim=2;
      size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
      unsigned long int ngridx=4*16;
      unsigned long int ngridy=(Ntuples/ngridx)+1;
      size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_Hexciton, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_Hexciton failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
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
      PrintOpenCLerrorMessage(err);
      exit(1);
   }
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      PrintOpenCLerrorMessage(err);
      exit(1);
   }
   double run_time =(double)(end_time - start_time);
   fprintf(stderr,"kernel_execute_Hexciton: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


void Kernel_HexcitonLargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ntuples, int Nsites, myreal hbar, myreal dt, cl_myreal2* Ham, OpenCL_Init* OpenCLinfo)
{
   int err;
   //
   cl_event event;
   if(strcmp(OpenCLinfo->kernel_type,"LoopIn")==0)
   {
      int workdim=1;
      size_t globalWorkSize[]= {Ntuples};
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,0,sizeof(INsigmanew), (void*)&INsigmanew);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,1,sizeof(INsigmaold), (void*)&INsigmaold);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,2,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,3,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,4,sizeof(myreal), (void*)&hbar);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,5,sizeof(myreal), (void*)&dt);
      clSetKernelArg(OpenCLinfo->kernel_execute_Hexciton,6,sizeof(Ham), (void*)&Ham);
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_Hexciton, workdim, 0, globalWorkSize, NULL,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_Hexciton failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
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
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,0,sizeof(INsigmanew), (void*)&INsigmanew);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,1,sizeof(INsigmaold), (void*)&INsigmaold);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,2,sizeof(int), (void*)&Ntuples);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,3,sizeof(int), (void*)&Nsites);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,4,sizeof(myreal), (void*)&hbar);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,5,sizeof(myreal), (void*)&dt);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,6,sizeof(Ham), (void*)&Ham);
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,7,sizeof(int), (void*)&Nstrides);
      // Dynamically allocate local memory (allocated per workgroup)
      clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonLargeSystem,8,2*dimBlockx*dimBlocky*2*sizeof(myreal),NULL);   // local/shared Memory
      int workdim=2;
      size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
      unsigned long int ngridx=4*16;
      unsigned long int ngridy=(Ntuples/ngridx)+1;
      size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
      err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HexcitonLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
      if(err<0)
      {
         cout<<"Error: Execute kernel_execute_HexcitonLargeSystem failed"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
         exit(1);
      }
      err=clWaitForEvents(1, &event);
      if(err<0)
      {
         cout<<"Error: clWaitForEvents"<<endl;
         cout<<"Error-Code: "<<err<<endl;
         PrintOpenCLerrorMessage(err);
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
      PrintOpenCLerrorMessage(err);
      exit(1);
   }
   err=clGetEventProfilingInfo(event,CL_PROFILING_COMMAND_END,sizeof(cl_ulong),&end_time,&return_bytes);
   if(err<0)
   {
      cout<<"Error: clGetEventProfilingInfo (end)"<<endl;
      cout<<"Error-Code: "<<err<<endl;
      PrintOpenCLerrorMessage(err);
      exit(1);
   }
   double run_time =(double)(end_time - start_time);
   fprintf(stderr,"kernel_execute_HexcitonLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}



#endif

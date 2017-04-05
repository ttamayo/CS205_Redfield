#ifndef Kernel_HEXCITONREDFIELD_H
#define Kernel_HEXCITONREDFIELD_H


// some Notes
//
// BEACHTE KOMPLEXE MULTIPLIKATIONSARITHMETIC RE=RE*RE-IM*IM, IM=RE*IM+IM*RE !!!!!!!
// Hexciton dynamic getestet


void Kernel_HexcitonRedfield(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Nsites, myreal hbar, myreal dt, cl_myreal2* Ham, OpenCL_Init* OpenCLinfo)
{
// definiere groesse der Thred-Bloecke
//  int ydim=((Nsites*Nsites)/16)+1*(Nsites*Nsites)%16; //define threadblock with dimension 16 x ydim, try to get muliples of 16, this increases performance of GPU
//  int nblockx=Nsites;//+Nsites%4;
//  int nblocky=Nsites;//+Nsites%2;
   int Ntuples=1; //in Redfield per default Ntuples=1
   int err;
   unsigned long int N=1*Nsites*Nsites;

   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,6,sizeof(Ham), (void*)&Ham);
   //
   cl_event event;

   // Dynamically allocate local memory (allocated per workgroup)
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfield,7,2*Nsites*Nsites*2*sizeof(myreal),NULL);   // local/shared Memory
   unsigned long int dimBlockx=Nsites;
   unsigned long int dimBlocky=Nsites;
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   // for Refield only one density matrix
   unsigned long int ngridx=1;
   unsigned long int ngridy=1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HexcitonRedfield, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_HexcitonRedfield failed"<<endl;
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
   fprintf(stderr,"kernel_execute_HexcitonRedfield: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


void Kernel_HexcitonRedfieldLargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Nsites, myreal hbar, myreal dt, cl_myreal2* Ham, OpenCL_Init* OpenCLinfo)
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
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,0,sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,1,sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,2,sizeof(int), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,3,sizeof(int), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,4,sizeof(myreal), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,5,sizeof(myreal), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,6,sizeof(Ham), (void*)&Ham);
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,7,sizeof(int), (void*)&Nstrides);
   // Dynamically allocate local memory (allocated per workgroup)
   clSetKernelArg(OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem,8,2*dimBlockx*dimBlocky*2*sizeof(myreal),NULL);   // local/shared Memory
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=1;
   unsigned long int ngridy=1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HexcitonRedfieldLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: Execute kernel_execute_HexcitonRedfieldLargeSystem failed"<<endl;
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
   fprintf(stderr,"kernel_execute_HexcitonRedfieldLargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}



#endif

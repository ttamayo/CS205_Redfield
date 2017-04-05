#ifndef GPU_HPhononDoubleExcitonLowTemp_H
#define GPU_HPhononDoubleExcitonLowTemp_H

void Kernel_HPhononDoubleLowTemp_Part1(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int* AllSigmaTuples, int* d_MemIDnminus1, int* d_MemIDnplus1, int* d_Memnplus1start, int* d_TupleNorm, int Nsites, int NsitesCoupled, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar ,myreal *d_S0, myreal *d_Imchi0, myreal *d_Rechi0, cl_myreal2* d_prefactLowTemp, cl_myreal2* d_prefactLowTemp2,  int* d_listNDL,  int* d_listNDLsite, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 2, sizeof(AllSigmaTuples), (void*)&AllSigmaTuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 3, sizeof(d_MemIDnminus1), (void*)&d_MemIDnminus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 4, sizeof(d_MemIDnplus1), (void*)&d_MemIDnplus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 5, sizeof(d_Memnplus1start), (void*)&d_Memnplus1start);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 6, sizeof(d_TupleNorm), (void*)&d_TupleNorm);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 7, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 8, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 9, sizeof(Ntuples), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 10, sizeof(Nmax), (void*)&Nmax);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 11, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 12, sizeof(hbar), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 13, sizeof(d_S0), (void*)&d_S0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 14, sizeof(d_Imchi0), (void*)&d_Imchi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 15, sizeof(d_Rechi0), (void*)&d_Rechi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 16, sizeof(d_prefactLowTemp), (void*)&d_prefactLowTemp);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 17, sizeof(d_prefactLowTemp2), (void*)&d_prefactLowTemp2);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 18, sizeof(d_listNDL), (void*)&d_listNDL);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, 19, sizeof(d_listNDLsite), (void*)&d_listNDLsite);
   //
   cl_event event;
   unsigned long int dimBlockx=NsitesCoupled;
   unsigned long int dimBlocky=Nsites;
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=4*16;
   unsigned int ngridy=(Ntuples/ngridx)+1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDoubleLowTemp_Part1"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDoubleLowTemp_Part1: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


void Kernel_HPhononDoubleLowTemp_Part1LargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold,  int* AllSigmaTuples,  int* d_MemIDnminus1,  int* d_MemIDnplus1,  int* d_Memnplus1start,  int* d_TupleNorm, int Nsites, int NsitesCoupled, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar ,myreal *d_S0, myreal *d_Imchi0, myreal *d_Rechi0, cl_myreal2* d_prefactLowTemp, cl_myreal2* d_prefactLowTemp2,  int* d_listNDL,  int* d_listNDLsite, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 2, sizeof(AllSigmaTuples), (void*)&AllSigmaTuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 3, sizeof(d_MemIDnminus1), (void*)&d_MemIDnminus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 4, sizeof(d_MemIDnplus1), (void*)&d_MemIDnplus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 5, sizeof(d_Memnplus1start), (void*)&d_Memnplus1start);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 6, sizeof(d_TupleNorm), (void*)&d_TupleNorm);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 7, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 8, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 9, sizeof(Ntuples), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 10, sizeof(Nmax), (void*)&Nmax);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 11, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 12, sizeof(hbar), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 13, sizeof(d_S0), (void*)&d_S0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 14, sizeof(d_Imchi0), (void*)&d_Imchi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 15, sizeof(d_Rechi0), (void*)&d_Rechi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 16, sizeof(d_prefactLowTemp), (void*)&d_prefactLowTemp);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 17, sizeof(d_prefactLowTemp2), (void*)&d_prefactLowTemp2);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 18, sizeof(d_listNDL), (void*)&d_listNDL);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 19, sizeof(d_listNDLsite), (void*)&d_listNDLsite);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, 20, sizeof(NsitesCoupled), (void*)&NsitesCoupled);
   cl_event event;
   unsigned long int dimBlockx=Nsites;
   int workdim=1;
   size_t localWorkSize[]= {dimBlockx}; // defines threadblock
   unsigned long int ngridx=Ntuples;
   size_t globalWorkSize[]= {ngridx*dimBlockx};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}



void Kernel_HPhononDouble_DoubleExcitationLowTemp(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int* AllSigmaTuples,  int* d_MemIDnminus1,  int* d_MemIDnplus1,  int* d_Memnplus1start,  int* d_TupleNorm, int Nsites, int Nsites1, int NsitesCoupled, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, myreal *d_S0, myreal *d_Imchi0, myreal *d_Rechi0, cl_myreal2* d_prefactLowTemp, cl_myreal2* d_prefactLowTemp2,  int* d_listNDL,  int* d_listNDLsite, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 2, sizeof(AllSigmaTuples), (void*)&AllSigmaTuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 3, sizeof(d_MemIDnminus1), (void*)&d_MemIDnminus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 4, sizeof(d_MemIDnplus1), (void*)&d_MemIDnplus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 5, sizeof(d_Memnplus1start), (void*)&d_Memnplus1start);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 6, sizeof(d_TupleNorm), (void*)&d_TupleNorm);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 7, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 8, sizeof(Nsites1), (void*)&Nsites1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 9, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 10, sizeof(Ntuples), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 11, sizeof(Nmax), (void*)&Nmax);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 12, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 13, sizeof(hbar), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 14, sizeof(d_S0), (void*)&d_S0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 15, sizeof(d_Imchi0), (void*)&d_Imchi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 16, sizeof(d_Rechi0), (void*)&d_Rechi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 17, sizeof(d_prefactLowTemp), (void*)&d_prefactLowTemp);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 18, sizeof(d_prefactLowTemp2), (void*)&d_prefactLowTemp2);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 19, sizeof(d_listNDL), (void*)&d_listNDL);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, 20, sizeof(d_listNDLsite), (void*)&d_listNDLsite);
   //
   cl_event event;
   unsigned long int dimBlockx=NsitesCoupled;
   unsigned long int dimBlocky=Nsites;
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=4*16;
   unsigned int ngridy=(Ntuples/ngridx)+1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTemp, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDouble_DoubleExcitationLowTemp"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDouble_DoubleExcitationLowTemp: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}
//
//
void Kernel_HPhononDouble_DoubleExcitationLowTempLargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int* AllSigmaTuples,  int* d_MemIDnminus1,  int* d_MemIDnplus1,  int* d_Memnplus1start,  int* d_TupleNorm, int Nsites, int Nsites1, int NsitesCoupled, int Ncoupled, int Ntuples, int Nmax, myreal dt, myreal hbar, myreal *d_S0, myreal *d_Imchi0, myreal *d_Rechi0, cl_myreal2* d_prefactLowTemp, cl_myreal2* d_prefactLowTemp2,  int* d_listNDL,  int* d_listNDLsite, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 2, sizeof(AllSigmaTuples), (void*)&AllSigmaTuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 3, sizeof(d_MemIDnminus1), (void*)&d_MemIDnminus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 4, sizeof(d_MemIDnplus1), (void*)&d_MemIDnplus1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 5, sizeof(d_Memnplus1start), (void*)&d_Memnplus1start);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 6, sizeof(d_TupleNorm), (void*)&d_TupleNorm);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 7, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 8, sizeof(Nsites1), (void*)&Nsites1);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 9, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 10, sizeof(Ntuples), (void*)&Ntuples);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 11, sizeof(Nmax), (void*)&Nmax);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 12, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 13, sizeof(hbar), (void*)&hbar);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 14, sizeof(d_S0), (void*)&d_S0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 15, sizeof(d_Imchi0), (void*)&d_Imchi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 16, sizeof(d_Rechi0), (void*)&d_Rechi0);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 17, sizeof(d_prefactLowTemp), (void*)&d_prefactLowTemp);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 18, sizeof(d_prefactLowTemp2), (void*)&d_prefactLowTemp2);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 19, sizeof(d_listNDL), (void*)&d_listNDL);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 20, sizeof(d_listNDLsite), (void*)&d_listNDLsite);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, 21, sizeof(NsitesCoupled), (void*)&NsitesCoupled);
   //
   cl_event event;
   unsigned long int dimBlockx=Nsites;
   int workdim=1;
   size_t localWorkSize[]= {dimBlockx}; // defines threadblock
   unsigned long int ngridx=Ntuples;
   size_t globalWorkSize[]= {ngridx*dimBlockx};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDouble_DoubleExcitationLowTemp"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDouble_DoubleExcitationLowTemp: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}




//
//
void Kernel_HPhononDoubleLowTemp_Part2(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ncoupled, int Nsites, int Ntuples, myreal dt, cl_myreal2* d_Tupleprefact, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 2, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 3, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 4, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 5, sizeof(d_Tupleprefact), (void*)&d_Tupleprefact);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, 6, sizeof(Ntuples), (void*)&Ntuples);
   cl_event event;
   unsigned long int dimBlockx=Nsites;
   unsigned long int dimBlocky=Nsites;
   int workdim=2;
   size_t localWorkSize[]= {dimBlockx,dimBlocky}; // defines threadblock
   unsigned long int ngridx=4*16;
   unsigned int ngridy=(Ntuples/ngridx)+1;
   size_t globalWorkSize[]= {ngridx*dimBlockx,ngridy*dimBlocky};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDouble_DoubleExcitationLowTemp"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDouble_DoubleExcitationLowTemp: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}


//
void Kernel_HPhononDoubleLowTemp_Part2LargeSystem(cl_myreal2* INsigmanew, cl_myreal2* INsigmaold, int Ncoupled, int Nsites, int Ntuples, myreal dt, cl_myreal2* d_Tupleprefact, OpenCL_Init* OpenCLinfo)
{
   int err;
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 0, sizeof(INsigmanew), (void*)&INsigmanew);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 1, sizeof(INsigmaold), (void*)&INsigmaold);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 2, sizeof(Ncoupled), (void*)&Ncoupled);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 3, sizeof(Nsites), (void*)&Nsites);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 4, sizeof(dt), (void*)&dt);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 5, sizeof(d_Tupleprefact), (void*)&d_Tupleprefact);
   clSetKernelArg(OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, 6, sizeof(Ntuples), (void*)&Ntuples);
   cl_event event;
   unsigned long int dimBlockx=Nsites;
   int workdim=1;
   size_t localWorkSize[]= {dimBlockx}; // defines threadblock
   unsigned long int ngridx=Ntuples;
   size_t globalWorkSize[]= {ngridx*dimBlockx};
   err=clEnqueueNDRangeKernel(OpenCLinfo->queue, OpenCLinfo->kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem, workdim, 0, globalWorkSize, localWorkSize,0,0,&event);
   if(err<0)
   {
      cout<<"Error: kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem"<<endl;
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
   fprintf(stderr,"kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem: %f ms\n",run_time/1.0e9*1000.0);
#endif
   clReleaseEvent(event);
}

#endif


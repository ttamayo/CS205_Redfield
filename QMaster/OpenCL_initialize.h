#ifndef OPENCL_INITIALIZE_H
#define OPENCL_INITIALIZE_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string.h>
#include <sstream>
#include <vector>
#include <cmath>
#include "headers.h"

using namespace std;

class OpenCL_Init
{
public:
   char* kernel_type;
   cl_platform_id *platforms;
   cl_uint num_platforms;
   cl_uint use_platform_vendor; // id of the platform of the specified vendor
   cl_device_id *devices;
   cl_uint use_device_id;
   cl_device_id *selected_device;
   cl_context context;
   cl_command_queue queue;  // context command queue
   cl_program program;
   cl_kernel kernel_execute_InitializeSigma;
   cl_kernel kernel_execute_aXsig1Plussig2;
   cl_kernel kernel_execute_ReturnNewSigma;
   cl_kernel kernel_execute_Hexciton;
   cl_kernel kernel_execute_HexcitonLargeSystem;
   cl_kernel kernel_execute_HexcitonRedfield;
   cl_kernel kernel_execute_HexcitonRedfieldLargeSystem;
   cl_kernel kernel_execute_antiHermitian;
   cl_kernel kernel_execute_antiHermitianLargeSystem;
   cl_kernel kernel_execute_antiHermitianRedfield;
   cl_kernel kernel_execute_antiHermitianRedfieldLargeSystem;
   cl_kernel kernel_execute_HPhononLowTemp_Part1;
   cl_kernel kernel_execute_HPhononLowTemp_Part2;
   cl_kernel kernel_execute_HPhononLowTemp_Part1LargeSystem;
   cl_kernel kernel_execute_HPhononLowTemp_Part2LargeSystem;
   cl_kernel kernel_execute_HPhononDoubleLowTemp_Part1;
   cl_kernel kernel_execute_HPhononDouble_DoubleExcitationLowTemp;
   cl_kernel kernel_execute_HPhononDoubleLowTemp_Part2;
   cl_kernel kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem;
   cl_kernel kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem;
   cl_kernel kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem;
   cl_kernel kernel_execute_SigmaMultiplicationLeft;
   cl_kernel kernel_execute_SigmaMultiplicationRight;
   cl_kernel kernel_execute_SigmaMultiplicationLeftLargeSystem;
   cl_kernel kernel_execute_SigmaMultiplicationRightLargeSystem;
   cl_kernel kernel_execute_sigmahelp_to_sigma;

   FILE *program_handle;
   char *program_log;
   Kernel_Source Kernels;
public:
   OpenCL_Init();
   void clean_up(void);
   void get_platforms(void);
   void get_kernels(char* INkernel_type);
   void initialize(char* vendor, char* device_type, int INuse_device_id);
};

void OpenCL_Init::get_kernels(char* INkernel_type)
{
   kernel_type=INkernel_type;
   cout << "# Initialize Kernels ...   "<<endl;
#ifdef DEVELOPMENT
   // For development use kernels in Kernels/*.cl. else: use Hardcoded kernels in KernelSource.h
   const char *cl_Filenames[]=
   {
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/InitSigma.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/RungeKutta.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/Hexciton.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/H_antiHermitian.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/HPhononLowTemp.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/HPhononDoubleExcitonLowTemp.cl",
      "/n/home00/ckreisbeck/QMaster_v2/src/Kernels/multiply_dipole_sigmas.cl"
   };
#endif
   // Now Kernels are Hardcoded in KernelSource.h, Translate Kernel files *.cl to hardcoded character is done with hardcode_kernels.py
   const char options[] = "-DDOUBLE_PRECISION -cl-finite-math-only -cl-no-signed-zeros";
   int err;
   int numFiles=7;
#ifdef DEVELOPMENT
   char* program_buffer[numFiles];
   size_t program_size[numFiles];
#endif
   size_t log_size;
#ifdef DEVELOPMENT
   for(int i=0; i<numFiles; i++)
   {
      cout<<"# load Kernel-File:  "<<cl_Filenames[i]<<endl;
      program_handle = fopen(cl_Filenames[i], "r");
      if(program_handle == NULL)
      {
         cout<<"Error:  Couldn't find the program file"<<endl;
         exit(1);
      }
      fseek(program_handle, 0, SEEK_END);
      program_size[i] = ftell(program_handle);
      rewind(program_handle);
      program_buffer[i] = (char*)malloc(program_size[i]+1);
      program_buffer[i][program_size[i]] = '\0';
      fread(program_buffer[i], sizeof(char), program_size[i], program_handle);
      fclose(program_handle);
   }
   program = clCreateProgramWithSource(context, numFiles,(const char**)program_buffer, program_size, &err);
#else
   program = clCreateProgramWithSource(context, Kernels.numFiles, (const char**)Kernels.program_buffer, Kernels.program_size, &err);
#endif
   if(err<0)
   {
      fprintf(stderr,"clCreateProgramWithSource ERROR %d\n",err);
   }
   else
   {
      fprintf(stdout,"# clCreateProgramWithSource SUCCESS\n");
   }
   err = clBuildProgram(program, 1, selected_device, options, NULL, NULL);
   if(err<0)
   {
      fprintf(stderr,"clBuildProgram ERROR %d\n",err);
      size_t log_size;
      clGetProgramBuildInfo(program, *selected_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
      // Allocate memory for the log
      char *log = (char *) malloc(log_size);
      // Get the log
      clGetProgramBuildInfo(program, *selected_device, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
      // Print the log
      fprintf(stderr,"%s\n", log);
      exit(1);
   }
   else
   {
      fprintf(stdout,"# clBuildProgram SUCCESS\n");
   }

   if(err < 0)
   {
      clGetProgramBuildInfo(program, *selected_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
      program_log = (char*)malloc(log_size+1);
      program_log[log_size] = '\0';
      clGetProgramBuildInfo(program, *selected_device,CL_PROGRAM_BUILD_LOG,log_size, program_log, NULL);
      printf("%s\n", program_log);
      free(program_log);
      exit(1);
   }
#ifdef DEVELOPMENT
   for(int i=0; i<numFiles; i++)
   {
      free(program_buffer[i]);
   }
#endif
/// Create Kernel
   // SIGMA_INIT
//   cout<<"# build SIGMA_INIT"<<endl;
   kernel_execute_InitializeSigma= clCreateKernel(program, "execute_InitializeSigma", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_InitializeSigma ERROR %d\n",err);
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_InitializeSigma SUCCESS\n");
   }
   // RUNGE_KUTTA
//   cout<<"# build RUNGE_KUTTA"<<endl;
   kernel_execute_aXsig1Plussig2= clCreateKernel(program, "execute_aXsig1Plussig2", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_aXsig1Plussig2 ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_aXsig1Plussig2 SUCCESS\n");
   }
   kernel_execute_ReturnNewSigma= clCreateKernel(program, "execute_ReturnNewSigma", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_ReturnNewSigma ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_ReturnNewSigma SUCCESS\n");
   }
   // Liouville Hexciton
//   cout<<"# build HEXCITON"<<endl;
   if(strcmp(kernel_type,"default")==0)
   {
      kernel_execute_Hexciton= clCreateKernel(program, "execute_Hexciton", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_Hexciton ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_Hexciton SUCCESS\n");
      }
      kernel_execute_HexcitonLargeSystem= clCreateKernel(program, "execute_HexcitonLargeSystem", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_HexcitonLargeSystem ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_HexcitonLargeSystem SUCCESS\n");
      }
   }
   if(strcmp(kernel_type,"LoopIn")==0)
   {
      kernel_execute_Hexciton= clCreateKernel(program, "execute_Hexciton_loopIn", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_Hexciton ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_Hexciton SUCCESS\n");
      }
   }
   // Liouville HexcitonRedfield
//   cout<<"# build HEXCITONREDFIELD"<<endl;
   kernel_execute_HexcitonRedfield= clCreateKernel(program, "execute_Hexciton", &err); //Refield use default Hexciton Kernel
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HexcitonRedfield ERROR %d\n",err);
      exit(1);
   }
   else
   {
//         fprintf(stdout,"# clCreateKernel kernel_execute_HexcitonRedfield SUCCESS\n");
   }
   kernel_execute_HexcitonRedfieldLargeSystem= clCreateKernel(program, "execute_HexcitonLargeSystem", &err); //Refield use default Hexciton Kernel
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HexcitonLargeSystem ERROR %d\n",err);
      exit(1);
   }
   else
   {
//         fprintf(stdout,"# clCreateKernel kernel_execute_HexcitonRedfieldLargeSystem SUCCESS\n");
   }

   // Liouville H_antiHermitian
//   cout<<"# build H_ANTIHERMITIAN"<<endl;
   if(strcmp(kernel_type,"default")==0)
   {
      kernel_execute_antiHermitian= clCreateKernel(program, "execute_antiHermitian", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_antiHermitian ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_antiHermitian SUCCESS\n");
      }
      kernel_execute_antiHermitianLargeSystem= clCreateKernel(program, "execute_antiHermitianLargeSystem", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_antiHermitianLargeSystem ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_antiHermitianLargeSystem SUCCESS\n");
      }
   }
   if(strcmp(kernel_type,"LoopIn")==0)
   {
      kernel_execute_antiHermitian= clCreateKernel(program, "execute_antiHermitian_loopIn", &err);
      if(err<0)
      {
         fprintf(stderr,"clCreateKernel kernel_execute_antiHermitian ERROR %d\n",err);
         exit(1);
      }
      else
      {
//         fprintf(stdout,"# clCreateKernel kernel_execute_antiHermitian SUCCESS\n");
      }
   }

   // Liouville H_antiHermitianRedfield
//   cout<<"# build H_ANTIHERMITIANREDFIELD"<<endl;
   kernel_execute_antiHermitianRedfield= clCreateKernel(program, "execute_antiHermitian", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_antiHermitianRedfield ERROR %d\n",err);
      exit(1);
   }
   else
   {
//         fprintf(stdout,"# clCreateKernel kernel_execute_antiHermitianRedfield SUCCESS\n");
   }
   kernel_execute_antiHermitianRedfieldLargeSystem= clCreateKernel(program, "execute_antiHermitianLargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_antiHermitianRedfieldLargeSystem ERROR %d\n",err);
      exit(1);
   }
   else
   {
//         fprintf(stdout,"# clCreateKernel kernel_execute_antiHermitianRedfieldLargeSystem SUCCESS\n");
   }

   // LIOUVILLE PHONONLOWTEMP
//   cout<<"# build PHONON"<<endl;
   kernel_execute_HPhononLowTemp_Part1= clCreateKernel(program, "execute_HPhononLowTemp_Part1", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononLowTemp_Part1 ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononLowTemp_Part1 SUCCESS\n");
   }
   kernel_execute_HPhononLowTemp_Part2= clCreateKernel(program, "execute_HPhononLowTemp_Part2", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononLowTemp_Part2 ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononLowTemp_Part2 SUCCESS\n");
   }
   kernel_execute_HPhononLowTemp_Part1LargeSystem= clCreateKernel(program, "execute_HPhononLowTemp_Part1LargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhonon_Part1LargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononLowTemp_Part1LargeSystem SUCCESS\n");
   }
   kernel_execute_HPhononLowTemp_Part2LargeSystem= clCreateKernel(program, "execute_HPhononLowTemp_Part2LargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononLowTemp_Part2LargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononLowTemp_Part2LargeSystem SUCCESS\n");
   }
//   fprintf(stdout, "# ... done\n");
   // LIOUVILLE PHONONDOUBLEEXCITONLOWTEMP
//   cout<<"# build PHONON-DOUBLE-EXCITON"<<endl;
   kernel_execute_HPhononDoubleLowTemp_Part1= clCreateKernel(program, "execute_HPhononDoubleLowTemp_Part1", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part1 ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part1 SUCCESS\n");
   }
   kernel_execute_HPhononDoubleLowTemp_Part2= clCreateKernel(program, "execute_HPhononDoubleLowTemp_Part2", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part2 ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part2 SUCCESS\n");
   }
   kernel_execute_HPhononDouble_DoubleExcitationLowTemp = clCreateKernel(program, "execute_HPhononDouble_DoubleExcitationLowTemp", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDouble_DoubleExcitationLowTemp ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDouble_DoubleExcitationLowTemp SUCCESS\n");
   }
   kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem= clCreateKernel(program, "execute_HPhononDoubleLowTemp_Part1LargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part1LargeSystem SUCCESS\n");
   }
   kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem= clCreateKernel(program, "execute_HPhononDoubleLowTemp_Part2LargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDoubleLowTemp_Part2LargeSystem SUCCESS\n");
   }
   kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem = clCreateKernel(program, "execute_HPhononDouble_DoubleExcitationLowTempLargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_HPhononDouble_DoubleExcitationLowTempLargeSystem SUCCESS\n");
   }
   // multiply_dipole_sigmas
//   cout<<"# build MULTIPLY_DIPOLE_SIGMAS"<<endl;
   kernel_execute_SigmaMultiplicationLeft= clCreateKernel(program, "execute_SigmaMultiplicationLeft", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_SigmaMultiplicationLeft ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_SigmaMultiplicationLeft SUCCESS\n");
   }
   kernel_execute_SigmaMultiplicationRight= clCreateKernel(program, "execute_SigmaMultiplicationRight", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_SigmaMultiplicationRight ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_SigmaMultiplicationRight SUCCESS\n");
   }
   kernel_execute_SigmaMultiplicationLeftLargeSystem= clCreateKernel(program, "execute_SigmaMultiplicationLeftLargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_SigmaMultiplicationLeftLargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_SigmaMultiplicationLeftLargeSystem SUCCESS\n");
   }
   kernel_execute_SigmaMultiplicationRightLargeSystem= clCreateKernel(program, "execute_SigmaMultiplicationRightLargeSystem", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_SigmaMultiplicationRightLargeSystem ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_SigmaMultiplicationRight SUCCESS\n");
   }
   kernel_execute_sigmahelp_to_sigma= clCreateKernel(program, "execute_sigmahelp_to_sigma", &err);
   if(err<0)
   {
      fprintf(stderr,"clCreateKernel kernel_execute_sigmahelp_to_sigma ERROR %d %s\n",err,oclErrorString(err));
      exit(1);
   }
   else
   {
//      fprintf(stdout,"# clCreateKernel kernel_execute_sigmahelp_to_sigma SUCCESS\n");
   }
   fprintf(stdout,"# clCreateKernels SUCCESS\n");
   fprintf(stdout, "# ... done\n");
   fflush(stdout);
}


OpenCL_Init::OpenCL_Init()
{
   use_platform_vendor=-1;
}


void OpenCL_Init::get_platforms()
{
   fprintf(stdout, "# Get information about OpenCL platforms ... \n");
   cl_uint num_devices; // numer of available devices
   size_t ext_size;
   char* plat_info; // platform information
   cl_uint num_entries ;
   cl_uint err_num_devices;
   cl_int err = clGetPlatformIDs(0, NULL, &num_platforms); ///  Initialize openCL Platforms
   if (err < 0)
   {
      fprintf(stderr, "# clGetPlatformIDs ERROR %d\n", err);
      fprintf(stderr, "# Error while searching for OpenCL platforms\n", err);
      exit(1);
   }
   fprintf(stdout, "# Number of available OpenCL platforms: %d\n",
           num_platforms);
   platforms = (cl_platform_id*) malloc(
                  sizeof(cl_platform_id) * num_platforms);
   clGetPlatformIDs(num_platforms, platforms, NULL);

   // Go through platforms and show information
   for (int i = 0; i < num_platforms; i++)
   {
      fprintf(stdout, "# ==> information about PLATFORM %d \n", i);
      clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &ext_size);
      plat_info = new char[ext_size];
      clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, ext_size, plat_info,NULL);
      cout << "#     " << plat_info << endl;
      delete[] plat_info;
      clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, 0, NULL, &ext_size);
      plat_info = new char[ext_size];
      clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, ext_size, plat_info,NULL);
      // cout << "#     " << plat_info << endl;
      delete[] plat_info;
      clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, 0, NULL,&ext_size);
      plat_info = new char[ext_size];
      clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, ext_size,plat_info, NULL);
      // cout << "#     " << plat_info << endl;
      delete[] plat_info;

      cl_device_id * loc_devices;
      clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, NULL, NULL,&num_devices);
      printf("#     Number of available OpenCL devices: %d\n", num_devices);
      loc_devices = (cl_device_id*) malloc(sizeof(cl_device_id) * num_devices);
      clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, num_devices, loc_devices, NULL);
      fflush(stdout);
      for(int i=0; i<num_devices; i++)
      {
         fprintf(stdout,"#     ==> information about DEVICE %d\n",i);
         size_t text;
         cl_uint buf_uint;
         cl_ulong buf_ulong;
         char* dev_info; // device information
         size_t workitem_dims;
         size_t workitem_size[3];
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_NAME, NULL, NULL, &text);
         dev_info = new char[text];
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_NAME, text, dev_info, NULL);
         printf("#         DEVICE_NAME = %s\n", dev_info);
         delete [] dev_info;
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_VENDOR, NULL, NULL, &text);
         dev_info = new char[text];
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_VENDOR, text, dev_info, NULL);
         // printf("#         DEVICE_VENDOR = %s\n", dev_info);
         delete [] dev_info;
         clGetDeviceInfo(loc_devices[i], CL_DRIVER_VERSION, NULL, NULL, &text);
         dev_info = new char[text];
         clGetDeviceInfo(loc_devices[i], CL_DRIVER_VERSION, text, dev_info, NULL);
         // printf("#         DRIVER_VERSION = %s\n", dev_info);
         delete [] dev_info;
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL);
         printf("#         DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
         // printf("#         DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
         // printf("#         DEVICE_GLOBAL_MEM_SIZE = %llu kB\n", (unsigned long long)buf_ulong/1024L);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
         // printf("#         DEVICE_MAX_WORK_GROUP_SIZE = %llu\n", (unsigned long long)buf_ulong);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(workitem_dims), &workitem_dims, NULL);
         // printf("#         DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%u\n", (unsigned int) workitem_dims);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
         // printf("#         DEVICE_MAX_WORK_ITEM_SIZES:\t%u / %u / %u \n", (unsigned int)workitem_size[0], (unsigned int)workitem_size[1], (unsigned int)workitem_size[2]);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
         // printf("#         DEVICE_LOCAL_MEM_SIZE = %llu kB\n", (unsigned long long)buf_ulong/1024L);
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_EXTENSIONS, NULL, NULL, &text);
         dev_info = new char[text];
         clGetDeviceInfo(loc_devices[i], CL_DEVICE_EXTENSIONS, text, dev_info, NULL);
         //printf("#         DEVICE_EXTENSIONS = %s\n", dev_info);
         delete [] dev_info;
         fflush(stdout);
      }
   }
   fprintf(stdout, "# ... done \n#\n");
   fflush(stdout);
}

void OpenCL_Init::initialize(char* vendor, char* device_type, int INuse_device_id)
{
   int device_id=INuse_device_id;
   cl_device_type use_device_type;
   fprintf(stdout, "# Initialize OpenCL device  ... \n");
   for (int i = 0; i < num_platforms; i++)
   {
      size_t ext_size;
      char* plat_info; // platform information
      clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &ext_size);
      plat_info = new char[ext_size];
      clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, ext_size, plat_info, NULL);
      string s1=vendor;
      string s2=plat_info;
      if (s2.find(s1) != std::string::npos)
      {
         use_platform_vendor=i;
         cout<<"# select platform vendor: "<<plat_info<<endl;
      }
   }
   if(use_platform_vendor==-1)
   {
      fprintf(stderr,"no vendor=%s found\n", vendor);
      fprintf(stderr,"found vendors are:\n");
      for (int i = 0; i < num_platforms; i++)
      {
         size_t ext_size;
         char* plat_info; // platform information
         clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &ext_size);
         plat_info = new char[ext_size];
         clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, ext_size, plat_info, NULL);
         fprintf(stderr,"%s\n",plat_info);
      }
      exit(1);
   }
   if(strcmp(device_type,"CPU")==0)
   {
      use_device_type=CL_DEVICE_TYPE_CPU;
   }
   else if(strcmp(device_type,"GPU")==0)
   {
      use_device_type=CL_DEVICE_TYPE_GPU;
   }
   else if(strcmp(device_type,"ACCELERATOR")==0)
   {
      use_device_type=CL_DEVICE_TYPE_ACCELERATOR;
   }
   else
   {
      fprintf(stderr,"error unknown 'device_type \n",device_type);
      fprintf(stderr,"possible choices: device_type=CPU/ACCELERATOR/GPU. ABORT.\n");
      exit(1);
   }
   cl_uint num_devices;   // numer of available devices
   clGetDeviceIDs(platforms[use_platform_vendor], use_device_type, NULL, NULL, &num_devices);
   devices = (cl_device_id*) malloc(sizeof(cl_device_id) * num_devices);
   clGetDeviceIDs(platforms[use_platform_vendor], use_device_type, num_devices, devices, NULL);
   cout<<"# select device_type: "<<device_type<<endl;
   if(num_devices==0)
   {
      bool found_alternative=false;
      fprintf(stderr,"no device_type=%s found on platform vendor=%s \n",device_type, vendor);
      fprintf(stderr,"search all other vendors for device_type=%s ... \n",device_type);
      for(int i=0; i<num_platforms; i++)
      {
         cl_uint num_devs=0;
         clGetDeviceIDs(platforms[i], use_device_type, NULL, NULL, &num_devs);
         if(num_devs>0)
         {
            size_t ext_size;
            char* plat_info; // platform information
            clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &ext_size);
            plat_info = new char[ext_size];
            clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, ext_size, plat_info, NULL);
            fprintf(stderr,"found %d devices of device_type=%s on vendor %s \n",num_devs, device_type, plat_info);
            found_alternative=true;
         }
      }
      fprintf(stderr,"... done \n",device_type);
      if(found_alternative==false)
      {
         fprintf(stderr,"no device_type=%s found on any vendor \n",device_type);
      }
      else
      {
         fprintf(stderr,"recommendation: use available device_type=%s on vendors listed above \n",device_type);
      }
      exit(1);
   }


   //PICK specific DEVICE
   // check if more than one device. If yes then demand that device_id in input file needs to be specified
   if(num_devices>1 and  device_id<0)
   {
      fprintf(stderr,"# ERROR: on vendor %s are multiple devices of type %s. Not specified which one should be used! \n",vendor, device_type);
      fprintf(stderr,"# ERROR: on vendor %s found %d devices of type %s \n", vendor, num_devices, device_type);
      fprintf(stderr,"# Recommendation: specify device_id in input file.\n# Possible choices are: \n", vendor, num_devices, device_type);
      for(int i=0; i<num_devices; i++)
      {
         fprintf(stderr,"# device_id=%d \n", i);
      }
      exit(1);
   }
   // if only one device, no need to specify device_id. use the available device
   else if(device_id<0)
   {
      device_id=0;
      use_device_id=device_id;
   }
   // if divice_id specified, check if there is a device with the specified id
   else if(device_id>num_devices-1)
   {
      fprintf(stderr,"# ERROR: on vendor %s number of found %d devices of type %s is smaller than the given device_id in the input file.\n", vendor, num_devices, device_type);
      fprintf(stderr,"# Recommendation: reduce specified device_id=%d in input file.\n# Possible choices are: \n",device_id);
      for(int i=0; i<num_devices; i++)
      {
         fprintf(stderr,"# device_id=%d \n", i);
      }
      exit(1);
   }
   // if suitable device_id specified, then use the specified device
   else
   {
      use_device_id=device_id;
   }
   // Print information of used device
   fprintf(stdout,"# ==> information about the selected DEVICE \n");
   char* dev_info;
   cl_uint buf_uint;
   cl_ulong buf_ulong;
   size_t workitem_dims;
   size_t workitem_size[3];
   //
   size_t text;
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_NAME, NULL, NULL, &text);
   dev_info = new char[text];
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_NAME, text, dev_info, NULL);
   printf("#     DEVICE_NAME = %s\n", dev_info);
   delete [] dev_info;
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL);
   printf("#     DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL);
   printf("#     DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
   printf("#     DEVICE_GLOBAL_MEM_SIZE = %llu kB\n", (unsigned long long)buf_ulong/1024L);
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
   // printf("#     DEVICE_MAX_WORK_GROUP_SIZE = %llu\n", (unsigned long long)buf_ulong);
   //clGetDeviceInfo(devices[use_device_id], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(workitem_dims), &workitem_dims, NULL);
   // printf("#     DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%u\n", (unsigned int) workitem_dims);
   // clGetDeviceInfo(devices[use_device_id], CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
   // printf("#     DEVICE_MAX_WORK_ITEM_SIZES:\t%u / %u / %u \n", (unsigned int)workitem_size[0], (unsigned int)workitem_size[1], (unsigned int)workitem_size[2]);
   clGetDeviceInfo(devices[use_device_id], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL);
   printf("#     DEVICE_LOCAL_MEM_SIZE = %llu kB\n", (unsigned long long)buf_ulong/1024L);
   //clGetDeviceInfo(devices[use_device_id], CL_DEVICE_EXTENSIONS, sizeof(buffer), buffer, NULL);
   // printf("#         DEVICE_EXTENSIONS = %s\n", buffer);
   //
   cl_int err;
   context = clCreateContext(NULL, 1, &devices[use_device_id], NULL, NULL, &err);  //use_device_id
   if(err<0)
   {
      fprintf(stderr,"clCreateContext ERROR %d\n",err);
   }
   else
   {
      fprintf(stdout,"# clCreateContext SUCCESS\n");
   }
   // Create a command queue
   // Create a command queue
#ifdef PROFILING
   queue = clCreateCommandQueue(context, devices[use_device_id], CL_QUEUE_PROFILING_ENABLE, &err);
#else
   queue = clCreateCommandQueue(context, devices[use_device_id], 0, &err);
#endif
   if(err<0)
   {
      fprintf(stderr,"clCreateCommandQueue ERROR %d\n",err);
   }
   else
   {
      fprintf(stdout,"# clCreateCommandQueue SUCCESS\n");
   }
   selected_device=&devices[use_device_id];
   //
   fprintf(stdout, "# ... done\n#\n");
}

void OpenCL_Init::clean_up()
{
   cout<<"# clear OpenCL context"<<endl;
   cl_int err;
   err= clFlush(queue);
   err= clFinish(queue);
   err= clReleaseCommandQueue(queue);
   err= clReleaseContext(context);
   delete [] platforms;
   delete [] devices;
#ifndef DEVELOPMENT
//   delete [] Kernels.program_buffer;
//   delete [] Kernels.program_size;
#endif
}
#endif

/* opencl_error.h */
#ifndef OPENCL_ERROR_H
#define OPENCL_ERROR_H

#define CL_SUCCESS 0
#define CL_DEVICE_NOT_FOUND -1
#define CL_DEVICE_NOT_AVAILABLE -2
#define CL_COMPILER_NOT_AVAILABLE -3
#define CL_MEM_OBJECT_ALLOCATION_FAILURE -4
#define CL_OUT_OF_RESOURCES -5
#define CL_OUT_OF_HOST_MEMORY -6
#define CL_PROFILING_INFO_NOT_AVAILABLE -7
#define CL_MEM_COPY_OVERLAP -8
#define CL_IMAGE_FORMAT_MISMATCH -9
#define CL_IMAGE_FORMAT_NOT_SUPPORTED -10
#define CL_BUILD_PROGRAM_FAILURE -11
#define CL_MAP_FAILURE -12

#define CL_INVALID_VALUE -30
#define CL_INVALID_DEVICE_TYPE -31
#define CL_INVALID_PLATFORM -32
#define CL_INVALID_DEVICE -33
#define CL_INVALID_CONTEXT -34
#define CL_INVALID_QUEUE_PROPERTIES -35
#define CL_INVALID_COMMAND_QUEUE -36
#define CL_INVALID_HOST_PTR -37
#define CL_INVALID_MEM_OBJECT -38
#define CL_INVALID_IMAGE_FORMAT_DESCRIPTOR -39
#define CL_INVALID_IMAGE_SIZE -40
#define CL_INVALID_SAMPLER -41
#define CL_INVALID_BINARY -42
#define CL_INVALID_BUILD_OPTIONS -43
#define CL_INVALID_PROGRAM -44
#define CL_INVALID_PROGRAM_EXECUTABLE -45
#define CL_INVALID_KERNEL_NAME -46
#define CL_INVALID_KERNEL_DEFINITION -47
#define CL_INVALID_KERNEL -48
#define CL_INVALID_ARG_INDEX -49
#define CL_INVALID_ARG_VALUE -50
#define CL_INVALID_ARG_SIZE -51
#define CL_INVALID_KERNEL_ARGS -52
#define CL_INVALID_WORK_DIMENSION -53
#define CL_INVALID_WORK_GROUP_SIZE -54
#define CL_INVALID_WORK_ITEM_SIZE -55
#define CL_INVALID_GLOBAL_OFFSET -56
#define CL_INVALID_EVENT_WAIT_LIST -57
#define CL_INVALID_EVENT -58
#define CL_INVALID_OPERATION -59
#define CL_INVALID_GL_OBJECT -60
#define CL_INVALID_BUFFER_SIZE -61
#define CL_INVALID_MIP_LEVEL -62
#define CL_INVALID_GLOBAL_WORK_SIZE -63
const char *PrintOpenCLerrorMessage(cl_int err)
{
   switch(err)
   {
   case 0:
      return "CL_SUCCESS";
   case -1:
      return "CL_DEVICE_NOT_FOUND";
   case -2:
      return "CL_DEVICE_NOT_AVAILABLE";
   case -3:
      return "CL_COMPILER_NOT_AVAILABLE";
   case -4:
      return "CL_MEM_OBJECT_ALLOCATION_FAILURE (out of device memory)";
   case -5:
      return "CL_OUT_OF_RESOURCES";
   case -6:
      return "CL_OUT_OF_HOST_MEMORY";
   case -7:
      return "CL_PROFILING_INFO_NOT_AVAILABLE";
   case -8:
      return "CL_MEM_COPY_OVERLAP";
   case -9:
      return "CL_IMAGE_FORMAT_MISMATCH";
   case -10:
      return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
   case -11:
      return "CL_BUILD_PROGRAM_FAILURE";
   case -12:
      return "CL_MAP_FAILURE";
   case -30:
      return "CL_INVALID_VALUE";
   case -31:
      return "CL_INVALID_DEVICE_TYPE";
   case -32:
      return "CL_INVALID_PLATFORM";
   case -33:
      return "CL_INVALID_DEVICE";
   case -34:
      return "CL_INVALID_CONTEXT";
   case -35:
      return "CL_INVALID_QUEUE_PROPERTIES";
   case -36:
      return "CL_INVALID_COMMAND_QUEUE";
   case -37:
      return "CL_INVALID_HOST_PTR";
   case -38:
      return "CL_INVALID_MEM_OBJECT";
   case -39:
      return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
   case -40:
      return "CL_INVALID_IMAGE_SIZE";
   case -41:
      return "CL_INVALID_SAMPLER";
   case -42:
      return "CL_INVALID_BINARY";
   case -43:
      return "CL_INVALID_BUILD_OPTIONS";
   case -44:
      return "CL_INVALID_PROGRAM";
   case -45:
      return "CL_INVALID_PROGRAM_EXECUTABLE";
   case -46:
      return "CL_INVALID_KERNEL_NAME";
   case -47:
      return "CL_INVALID_KERNEL_DEFINITION";
   case -48:
      return "CL_INVALID_KERNEL";
   case -49:
      return "CL_INVALID_ARG_INDEX";
   case -50:
      return "CL_INVALID_ARG_VALUE";
   case -51:
      return "CL_INVALID_ARG_SIZE";
   case -52:
      return "CL_INVALID_KERNEL_ARGS";
   case -53:
      return "CL_INVALID_WORK_DIMENSION";
   case -54:
      return "CL_INVALID_WORK_GROUP_SIZE";
   case -55:
      return "CL_INVALID_WORK_ITEM_SIZE";
   case -56:
      return "CL_INVALID_GLOBAL_OFFSET";
   case -57:
      return "CL_INVALID_EVENT_WAIT_LIST";
   case -58:
      return "CL_INVALID_EVENT";
   case -59:
      return "CL_INVALID_OPERATION";
   case -60:
      return "CL_INVALID_GL_OBJECT";
   case -61:
      return "CL_INVALID_BUFFER_SIZE";
   case -62:
      return "CL_INVALID_MIP_LEVEL";
   case -63:
      return "CL_INVALID_GLOBAL_WORK_SIZE";
   default:
      return "Unknown OpenCL error";
   }

}
//void PrintOpenCLerrorMessage(cl_int err)
//{
//        switch (err) {
//        case CL_SUCCESS:                            cout<< "Success!"<<endl;
//        case CL_DEVICE_NOT_FOUND:                   cout<< "Device not found."<<endl;
//        case CL_DEVICE_NOT_AVAILABLE:               cout<< "Device not available"<<endl;
//        case CL_COMPILER_NOT_AVAILABLE:             cout<< "Compiler not available"<<endl;
//        case CL_MEM_OBJECT_ALLOCATION_FAILURE:      cout<< "Memory object allocation failure"<<endl;
//        case CL_OUT_OF_RESOURCES:                   cout<< "Out of resources"<<endl;
//        case CL_OUT_OF_HOST_MEMORY:                 cout<< "Out of host memory"<<endl;
//        case CL_PROFILING_INFO_NOT_AVAILABLE:       cout<< "Profiling information not available"<<endl;
//        case CL_MEM_COPY_OVERLAP:                   cout<< "Memory copy overlap"<<endl;
//        case CL_IMAGE_FORMAT_MISMATCH:              cout<< "Image format mismatch"<<endl;
//        case CL_IMAGE_FORMAT_NOT_SUPPORTED:         cout<< "Image format not supported"<<endl;
//        case CL_BUILD_PROGRAM_FAILURE:              cout<< "Program build failure"<<endl;
//        case CL_MAP_FAILURE:                        cout<< "Map failure"<<endl;
//        case CL_INVALID_VALUE:                      cout<< "Invalid value"<<endl;
//        case CL_INVALID_DEVICE_TYPE:                cout<< "Invalid device type"<<endl;
//        case CL_INVALID_PLATFORM:                   cout<< "Invalid platform"<<endl;
//        case CL_INVALID_DEVICE:                     cout<< "Invalid device"<<endl;
//        case CL_INVALID_CONTEXT:                    cout<< "Invalid context"<<endl;
//        case CL_INVALID_QUEUE_PROPERTIES:           cout<< "Invalid queue properties"<<endl;
//        case CL_INVALID_COMMAND_QUEUE:              cout<< "Invalid command queue"<<endl;
//        case CL_INVALID_HOST_PTR:                   cout<< "Invalid host pointer"<<endl;
//        case CL_INVALID_MEM_OBJECT:                 cout<< "Invalid memory object"<<endl;
//        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    cout<< "Invalid image format descriptor"<<endl;
//        case CL_INVALID_IMAGE_SIZE:                 cout<< "Invalid image size"<<endl;
//        case CL_INVALID_SAMPLER:                    cout<< "Invalid sampler"<<endl;
//        case CL_INVALID_BINARY:                     cout<< "Invalid binary"<<endl;
//        case CL_INVALID_BUILD_OPTIONS:              cout<< "Invalid build options"<<endl;
//        case CL_INVALID_PROGRAM:                    cout<< "Invalid program"<<endl;
//        case CL_INVALID_PROGRAM_EXECUTABLE:         cout<< "Invalid program executable"<<endl;
//        case CL_INVALID_KERNEL_NAME:                cout<< "Invalid kernel name"<<endl;
//        case CL_INVALID_KERNEL_DEFINITION:          cout<< "Invalid kernel definition"<<endl;
//        case CL_INVALID_KERNEL:                     cout<< "Invalid kernel"<<endl;
//        case CL_INVALID_ARG_INDEX:                  cout<< "Invalid argument index"<<endl;
//        case CL_INVALID_ARG_VALUE:                  cout<< "Invalid argument value"<<endl;
//        case CL_INVALID_ARG_SIZE:                   cout<< "Invalid argument size"<<endl;
//        case CL_INVALID_KERNEL_ARGS:                cout<< "Invalid kernel arguments"<<endl;
//        case CL_INVALID_WORK_DIMENSION:             cout<< "Invalid work dimension"<<endl;
//        case CL_INVALID_WORK_GROUP_SIZE:            cout<< "Invalid work group size"<<endl;
//        case CL_INVALID_WORK_ITEM_SIZE:             cout<< "Invalid work item size"<<endl;
//        case CL_INVALID_GLOBAL_OFFSET:              cout<< "Invalid global offset"<<endl;
//        case CL_INVALID_EVENT_WAIT_LIST:            cout<< "Invalid event wait list"<<endl;
//        case CL_INVALID_EVENT:                      cout<< "Invalid event"<<endl;
//        case CL_INVALID_OPERATION:                  cout<< "Invalid operation"<<endl;
//        case CL_INVALID_GL_OBJECT:                  cout<< "Invalid OpenGL object"<<endl;
//        case CL_INVALID_BUFFER_SIZE:                cout<< "Invalid buffer size"<<endl;
//        case CL_INVALID_MIP_LEVEL:                  cout<< "Invalid mip-map level"<<endl;
//        default: cout<<"unspecified error, check OpenCL errors online"<<endl;
//    }
//        cout.flush();
//}



const char* oclErrorString(cl_int error)
{
   static const char* errorString[] =
   {
      "CL_SUCCESS",
      "CL_DEVICE_NOT_FOUND",
      "CL_DEVICE_NOT_AVAILABLE",
      "CL_COMPILER_NOT_AVAILABLE",
      "CL_MEM_OBJECT_ALLOCATION_FAILURE",
      "CL_OUT_OF_RESOURCES",
      "CL_OUT_OF_HOST_MEMORY",
      "CL_PROFILING_INFO_NOT_AVAILABLE",
      "CL_MEM_COPY_OVERLAP",
      "CL_IMAGE_FORMAT_MISMATCH",
      "CL_IMAGE_FORMAT_NOT_SUPPORTED",
      "CL_BUILD_PROGRAM_FAILURE",
      "CL_MAP_FAILURE",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "",
      "CL_INVALID_VALUE",
      "CL_INVALID_DEVICE_TYPE",
      "CL_INVALID_PLATFORM",
      "CL_INVALID_DEVICE",
      "CL_INVALID_CONTEXT",
      "CL_INVALID_QUEUE_PROPERTIES",
      "CL_INVALID_COMMAND_QUEUE",
      "CL_INVALID_HOST_PTR",
      "CL_INVALID_MEM_OBJECT",
      "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
      "CL_INVALID_IMAGE_SIZE",
      "CL_INVALID_SAMPLER",
      "CL_INVALID_BINARY",
      "CL_INVALID_BUILD_OPTIONS",
      "CL_INVALID_PROGRAM",
      "CL_INVALID_PROGRAM_EXECUTABLE",
      "CL_INVALID_KERNEL_NAME",
      "CL_INVALID_KERNEL_DEFINITION",
      "CL_INVALID_KERNEL",
      "CL_INVALID_ARG_INDEX",
      "CL_INVALID_ARG_VALUE",
      "CL_INVALID_ARG_SIZE",
      "CL_INVALID_KERNEL_ARGS",
      "CL_INVALID_WORK_DIMENSION",
      "CL_INVALID_WORK_GROUP_SIZE",
      "CL_INVALID_WORK_ITEM_SIZE",
      "CL_INVALID_GLOBAL_OFFSET",
      "CL_INVALID_EVENT_WAIT_LIST",
      "CL_INVALID_EVENT",
      "CL_INVALID_OPERATION",
      "CL_INVALID_GL_OBJECT",
      "CL_INVALID_BUFFER_SIZE",
      "CL_INVALID_MIP_LEVEL",
      "CL_INVALID_GLOBAL_WORK_SIZE",
   };
   const int errorCount = sizeof(errorString) / sizeof(errorString[0]);
   const int index = -error;
   return (index >= 0 && index < errorCount) ? errorString[index] : "Unspecified Error";
}


#endif


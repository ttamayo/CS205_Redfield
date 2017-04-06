/*
 * OpenCL program to add two square real matrices with float entries
 *     i.e.  A + B = C
 * 
 * author: hannahsim
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define MAX_SOURCE_SIZE (0x100000)

// Working with one device for now
#define MAX_PLATFORMS 1
#define MAX_DEVICES 1

// *******************************************************************
// PARAMETERS
// *******************************************************************
#define DATA_SIZE 32 // matrix size

#define KERNEL_FILE_NAME "matrix_addition.cl" // kernel file name
#define LOCAL_SIZE 16*16 // analogous to CUDA's thread group


// DON'T CHANGE
#define GLOBAL_SIZE DATA_SIZE*DATA_SIZE // analogous to CUDA's block group



// *******************************************************************
// FOR TESTING PURPOSES
// *******************************************************************

// Allocates a matrix with designated entries (1 for testing purposes).
void matrixInit(float* data, int ssize) {
    for (int i = 0; i < ssize; ++i)
    data[i] = 1;
}

// Function to check output matrix C (for testing)
// A and B each have all entries 1 so we know what C's entries should be.
float checkResult(float* data, int ssize) {
    float tot = 0;
    for (int i =0; i < ssize; ++i) {
        tot += data[i];
    }
    return tot;
}



// *******************************************************************
// DEFINE KERNEL
//
// Update: can either include kernel in this file or have
// 	   a separate file for the kernel. I decided to 
// 	   have another file: matrix_addition.cl
// *******************************************************************
/*
const char *source = {
    "__kernel void matrix_add(__global float* matAA, __global float* matBB, __global float* matCC)\n"
    "{\n"
    "    int i = get_global_id(0);\n"
    "    matCC[i] = matAA[i] + matBB[i];\n"
    "}\n"
};
*/




// *******************************************************************
// PROGRAM MAIN
// *******************************************************************
int main(void) {

    clock_t tic = clock();   

    // for processing file containing kernel
    FILE *fp;
    char *source;
    size_t source_size;
    fp = fopen(KERNEL_FILE_NAME, "r");
    if (!fp) {
        printf("Error: Failed to load kernel file.\n");
        return EXIT_FAILURE;
    }
    source = (char*)malloc(MAX_SOURCE_SIZE);
    source_size = fread( source, 1, MAX_SOURCE_SIZE, fp);
    fclose( fp );   
  
    unsigned int i; // iterator variable

    printf("\n\n*********************************************\n\n");
    printf("Matrix Addition Implementation using OpenCL:");
    printf("\n*********************************************\n");

    int size = DATA_SIZE;
    printf("Matrix size: %u\n\n",size);
    printf("Expected output (sum of result matrix entries): %f\n\n",2.*size*size); // testing


    // Allocate host memory for matrices A and B 
    unsigned int size_A = size * size;
    unsigned int mem_size_A = sizeof(float) * size_A;
    float* h_A = (float*) malloc(mem_size_A);

    unsigned int size_B = size * size;
    unsigned int mem_size_B = sizeof(float) * size_B;
    float* h_B = (float*) malloc(mem_size_B);
    //printf("Checkpoint: Allocated host memory for input matrices A and B. \n");


    // Initialize host memory
    matrixInit(h_A, size_A);
    matrixInit(h_B, size_B);
    //printf("Checkpoint: Initialized host memory for input matrices A and B. \n");


    // Allocate host memory for the result matrix C
    unsigned int size_C = size * size;
    unsigned int mem_size_C = sizeof(float) * size_C;
    float* h_C = (float*) malloc(mem_size_C);
    //printf("Checkpoint: Allocated host memory for output matrix C. \n");
    


    ///////////////////////////////////////////////
    //                Setup OpenCL               //
    ///////////////////////////////////////////////

    // OpenCL-specific variables
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;
    cl_int err;
    
    // OpenCL device memory (memory objects) for matrices
    cl_mem d_A; 
    cl_mem d_B;
    cl_mem d_C;

    // Variables for device(s) and platform(s)
    // For now, using 1 device and 1 platform
    cl_platform_id platforms[MAX_PLATFORMS]; // IDs of all platforms
    cl_uint num_platforms; // number of platforms on machine
    cl_device_id devices[MAX_DEVICES]; // IDs of devices for each platform
    cl_uint num_devices; // number of devices on machine


    // Get number of platforms
    err = clGetPlatformIDs(1, platforms, &num_platforms);
    if (num_platforms != 1) {
        printf("Number of platforms: %u\n",num_platforms);
    }

 
    // Get number of devices on particular platform
    // OPTION: Can use GPU specifically via replacing CL_DEVICE_TYPE_ALL to CL_DEVICE_TYPE_GPU
    for (i=0; i<num_platforms; i++) {
        err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, sizeof(devices), devices, &num_devices);
        if (err != CL_SUCCESS) exit((printf("Error in clGetDeviceIDs %d",err),1));
        if (num_devices != 1) {
            printf("Number of devices: %u\n",num_devices);
        }
    }


    // Create context (for sharing between devices)
    context = clCreateContext(NULL, num_devices, devices, NULL, NULL, &err);
    if (!context) {
        printf("Error: Failed to create a compute context ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Computed context created. \n");
 

    // Create a command queue for a device (to do work on device)
    // NOTE: A command queue per device (but we're only using 1 device)
    // IF we use more than 1 device, need to create multiple command queues...
    queue = clCreateCommandQueue(context, devices[0], 0, &err);
    if (!queue) {
        printf("Error: Failed to create a command queue ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: A command queue created for a device.\n");


    // Set up device memory (create memory objects)
    d_A = clCreateBuffer(context, CL_MEM_READ_ONLY, mem_size_A, NULL, &err);
    d_B = clCreateBuffer(context, CL_MEM_READ_ONLY, mem_size_B, NULL, &err);
    d_C = clCreateBuffer(context, CL_MEM_WRITE_ONLY, mem_size_C, NULL, &err);
    if (!d_A || !d_B || !d_C) {
        printf("Error: Failed to allocate device memory ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Created memory objects for matrices.\n");


    // Copy matrices A and B to respective memory buffers
    err = clEnqueueWriteBuffer(queue, d_A, CL_TRUE, 0, mem_size_A, h_A, 0, NULL, NULL);
    err |= clEnqueueWriteBuffer(queue, d_B, CL_TRUE, 0, mem_size_B, h_B, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to copy matrices A and B to memory buffers ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Copied input matrix data to memory buffers. \n");


    // Create OpenCL program with source code
    program = clCreateProgramWithSource(context, 1, 
            (const char **)&source, (const size_t *)&source_size, &err); // use this if separate file for kernel
    // program = clCreateProgramWithSource(context, 1, (const char**)&source, NULL, &err); // use this if kernel included in this code
    if (!program) {
        printf("Error: Failed to create compute program ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Compute program created.\n");


    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to build program executable ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Built program.\n");


    // Create kernel object in the program we wish to run
    kernel = clCreateKernel(program, "matrix_add", &err);
    if (!kernel || err != CL_SUCCESS) {
        printf("Error: Failed to create compute kernel ... \n");
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Created compute kernel.\n");

    
    // Set kernel arguments
    size_t local, global; 
  
    err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&d_A);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&d_B);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&d_C);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to set kernel arguments ... %d\n", err);
        return EXIT_FAILURE;
    }

    cl_uint work_dim = 1;
    local = LOCAL_SIZE;
    global = GLOBAL_SIZE;


    // Copy data to the device then execute kernel
    err = clEnqueueNDRangeKernel(queue, kernel, work_dim, NULL, &global, &local, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to execute kernel ...  %d\n", err);
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Executed kernel.\n");

 
    // Retrieve result from device
    err = clEnqueueReadBuffer(queue, d_C, CL_TRUE, 0, mem_size_C, h_C, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array ... %d\n", err);
        return EXIT_FAILURE;
    }
    //printf("Checkpoint: Retrieved result from device.\n");
   


/*
    // Display results
    // print out A and B 
    printf("\n\nMatrix A\n");
    for(int i = 0; i < size_A; i++) {
        printf("%f ", h_A[i]);
        if(((i + 1) % size) == 0)
        printf("\n");
    }

    printf("\n\nMatrix B\n");
    for(int i = 0; i < size_B; i++) {
        printf("%f ", h_B[i]);
        if(((i + 1) % size) == 0)
        printf("\n");
    }

    // Print out result matrix C
    printf("\n\nMatrix C (Results)\n");
    for(int i = 0; i < size_C; i++) {
        printf("%f ", h_C[i]);
        if(((i + 1) % size) == 0)
        printf("\n");
    }
    printf("\n");
*/


    // Check result (for testing purposes)
    float sumOfC = checkResult(h_C, size_C);
    printf("\nCHECK: Sum of entries of matrix C: %f\n",sumOfC);



    // Clean up
    clFlush(queue);
    clFinish(queue); // Wait for everything to finish...

    clReleaseMemObject(d_A);
    clReleaseMemObject(d_B);
    clReleaseMemObject(d_C);
    
    clReleaseContext(context);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
 
    free(h_A);
    free(h_B);
    free(h_C);

    free(source);

    // Time program.
    clock_t toc = clock();
    double time_spent = (double)(toc-tic) / CLOCKS_PER_SEC;
    printf("\nTime Elapsed: %f seconds\n\n", time_spent);
    printf("\n*********************************************\n\n"); 

    return 0;
 
}

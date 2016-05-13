#ifndef OPENCLENV_H
#define OPENCLENV_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <glog/logging.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

class OpenCLEnv 
: public boost::basic_lockable_adapter<boost::mutex> {
 public:
  OpenCLEnv(
      const char* bin_path,
      const char* kernel_name): initialized(false)
  {
    // start platform setting up
    int err;

    cl_platform_id platform_id;
    cl_device_id device_id;

    char cl_platform_vendor[1001];
    char cl_platform_name[1001];

    cl_uint num_platforms = 0;

    // Connect to first platform
    err = clGetPlatformIDs(2, &platform_id, &num_platforms);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "Failed to find an OpenCL platform!");
    }

    err = clGetPlatformInfo(
        platform_id, 
        CL_PLATFORM_VENDOR, 
        1000, 
        (void *)cl_platform_vendor,NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!");
    }
    err = clGetPlatformInfo(platform_id,CL_PLATFORM_NAME,1000,(void *)cl_platform_name,NULL);
    
    if (err != CL_SUCCESS) {
        throw std::runtime_error("clGetPlatformInfo(CL_PLATFORM_NAME) failed!");
    }

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group!");
    }

    // Create a compute context 
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);

    if (!context) {
        throw std::runtime_error("Failed to create a compute context!");
    }

    // Create a command commands
    commands = clCreateCommandQueue(context, device_id, 0, &err);

    if (!commands) {
        throw std::runtime_error("Failed to create a command queue context!");
    }

    // Create Program Objects
    // TODO: this part should not be static, the program
    // should be configurable at runtime

    // Load binary from disk
    unsigned char *kernelbinary;
    
    int n_i = load_file(bin_path, (char **) &kernelbinary);

    if (n_i < 0) {
      throw std::runtime_error(
          "failed to load kernel from xclbin");
    }
    size_t n_t = n_i;

    int status = 0;

    // Create the compute program from offline
    program = clCreateProgramWithBinary(context, 1, &device_id, &n_t,
            (const unsigned char **) &kernelbinary, &status, &err);

    if ((!program) || (err!=CL_SUCCESS)) {
        throw std::runtime_error(
            "Failed to create compute program from binary");
    }

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to build program executable!");
    }

    // Create the compute kernel in the program we wish to run
    kernel = clCreateKernel(program, kernel_name, &err);

    if (!kernel || err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create compute kernel!");
    }

    initialized = true;
  }

  ~OpenCLEnv() {
    if (initialized) {
      clReleaseKernel(kernel);
      clReleaseCommandQueue(commands);
      clReleaseProgram(program);
      clReleaseContext(context);
    }
  }

  cl_context& getContext() {
    if (initialized) {
      return context;
    }
    else {
      throw std::runtime_error("environment not setup");
    }
  }

  cl_command_queue& getCmdQueue() {
    if (initialized) {
      return commands;
    }
    else {
      throw std::runtime_error("environment not setup");
    }
  }

  cl_kernel& getKernel() {
    if (initialized) {
      return kernel;
    }
    else {
      throw std::runtime_error("environment not setup");
    }
  }

 private:
  int load_file(
      const char *filename, 
      char **result)
  { 
    int size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) 
    { 
      *result = NULL;
      return -1; // -1 means file opening fail 
    } 
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) 
    { 
      free(*result);
      return -2; // -2 means file reading fail 
    } 
    fclose(f);
    (*result)[size] = 0;
    return size;
  }

  bool initialized;

  cl_context context;                 // compute context
  cl_command_queue commands;          // compute command queue
  cl_program program;                 // compute program
  cl_kernel kernel;                   // compute kernel
};
#endif

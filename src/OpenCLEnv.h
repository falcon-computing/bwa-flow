#ifndef OPENCLENV_H
#define OPENCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <glog/logging.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

#include "kflow/Queue.h"
#include "util.h"

class OpenCLEnv 
: public boost::basic_lockable_adapter<boost::mutex> {
 public:
  OpenCLEnv(const char* bin_path, const char* kernel_name) {
    // start platform setting up
    int err;

    cl_platform_id* platform_id = new cl_platform_id[2];

    char cl_platform_vendor[1001];
    char cl_platform_name[1001];

    cl_uint num_platforms = 0;

    // Connect to first platform
    err = clGetPlatformIDs(2, platform_id, &num_platforms);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "Failed to find an OpenCL platform!");
    }
    DLOG(INFO) << "Found " << num_platforms << " opencl platforms";

    int platform_idx;
    for (int i = 0; i < num_platforms; i++) {
        char cl_platform_name[1001];

        err = clGetPlatformInfo(
            platform_id[i], 
            CL_PLATFORM_NAME, 
            1000, 
            (void *)cl_platform_name,NULL);

        DLOG(INFO) << "Found platform " << cl_platform_name;

        if (err != CL_SUCCESS) {
          throw std::runtime_error(
              "clGetPlatformInfo(CL_PLATFORM_NAME) failed!");
        }

        if (strstr(cl_platform_name, "Xilinx") != NULL || 
            strstr(cl_platform_name, "Intel") != NULL ||
            strstr(cl_platform_name, "Altera") != NULL) {
            // found platform
            platform_idx = i;
            DLOG(INFO) << "Selecting platform " << cl_platform_name;
            break;
        }
    }

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id_, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group: " +
            std::to_string((long long)err));
    }

    // Create a compute context 
    context_ = clCreateContext(0, 1, &device_id_, NULL, NULL, &err);

    if (!context_) {
      throw std::runtime_error("Failed to create a compute context!");
    }

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
    program_ = clCreateProgramWithBinary(context_, 1, &device_id_, &n_t,
            (const unsigned char **) &kernelbinary, &status, &err);

    if ((!program_) || (err!=CL_SUCCESS)) {
        throw std::runtime_error(
            "Failed to create compute program from binary");
    }

    // Build the program executable
    err = clBuildProgram(program_, 0, NULL, NULL, NULL, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to build program executable!");
    }
  }

  ~OpenCLEnv() {
    clReleaseProgram(program_);
    clReleaseContext(context_);
  }

  cl_context& getContext() {
    return context_;
  }

  cl_device_id getDeviceId() {
    return device_id_;
  }

  cl_program getProgram() {
    return program_;
  }

 private:
  int load_file(
      const char *filename, 
      char **result) { 
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
  cl_context    context_;      // compute context
  cl_device_id  device_id_;
  cl_program    program_;      // compute program
};

extern OpenCLEnv* opencl_env;
#endif

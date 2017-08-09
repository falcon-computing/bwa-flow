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
#ifdef XILINX_FPGA
     // Create and write the reference buffer
     
     FILE* pac_inp = fopen("/pool/storage/yaoh/human_g1k_v37.fasta.pac", "rb");
     int pac_size = 0;
     fseek(pac_inp, 0, SEEK_END);
     pac_size = ftell(pac_inp);
     fseek(pac_inp, 0, SEEK_SET);
     pac_size = (pac_size+3)/4;

     int *pac = (int*)malloc(sizeof(int)*pac_size);
     fread(pac, sizeof(int), pac_size, pac_inp);
     fclose(pac_inp);

     command_ = clCreateCommandQueue(context_, device_id_, 0, &err);
     if (err != CL_SUCCESS) {
       throw std::runtime_error("Failed to create a command queue context!");
     }

     cl_mem_ext_ptr_t ext_c, ext_d;
     ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
     ext_d.flags = XCL_MEM_DDR_BANK3; ext_d.obj = 0; ext_d.param = 0;

     pac_input_a_ = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
         sizeof(int) * pac_size, &ext_c, &err);
     pac_input_b_ = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
         sizeof(int) * pac_size, &ext_d, &err);
     if (err != CL_SUCCESS) {
       throw std::runtime_error("Failed to create reference OpenCL buffer!");
     }
     cl_event event[2];
     err = clEnqueueWriteBuffer(command_, pac_input_a_, CL_TRUE, 0, sizeof(int) * pac_size, pac, 0, NULL, &event[0]);
     err = clEnqueueWriteBuffer(command_, pac_input_b_, CL_TRUE, 0, sizeof(int) * pac_size, pac, 0, NULL, &event[1]);
     clWaitForEvents(2, event);
     if (err != CL_SUCCESS) {
       throw std::runtime_error("Failed to write reference to DDR!");
     }
     free(pac);
     clReleaseEvent(event[0]);
     clReleaseEvent(event[1]);
#endif
   }

  ~OpenCLEnv() {
#ifdef XILINX_FPGA 
    clReleaseCommandQueue(command_);
    clReleaseMemObject(pac_input_a_);
    clReleaseMemObject(pac_input_b_);
#endif
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
#ifdef XILINX_FPGA
  cl_command_queue command_;
  cl_mem pac_input_a_;
  cl_mem pac_input_b_;
#endif
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

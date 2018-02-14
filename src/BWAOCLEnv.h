#ifndef BWAOCLENV_H
#define BWAOCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <glog/logging.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

#include "kflow/Queue.h"
#include "OpenCLEnv.h"
#include "util.h"

class BWAOCLEnv :public OpenCLEnv{
 public:
  BWAOCLEnv(const char* bin_path, const char* pac_path, const char* kernel_name)
    :OpenCLEnv(bin_path, kernel_name) {
     FILE* pac_inp = fopen(pac_path, "rb");
     int pac_size = 0;
     fseek(pac_inp, 0, SEEK_END);
     pac_size = ftell(pac_inp);
     fseek(pac_inp, 0, SEEK_SET);
     pac_size = (pac_size+3)/4;

     int *pac = (int*)malloc(sizeof(int)*pac_size);
     fread(pac, sizeof(int), pac_size, pac_inp);
     fclose(pac_inp);

     int err = 0;
     command_ = clCreateCommandQueue(context_, device_id_, 0, &err);
     if (err != CL_SUCCESS) {
       throw std::runtime_error("Failed to create a command queue context!");
     }
#ifdef XILINX_FPGA
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
#else
    DLOG(ERROR) << "This feature is currently only supported in Xilinx";
#endif
  }
  ~BWAOCLEnv(){
    clReleaseCommandQueue(command_);
#ifdef XILINX_FPGA
    clReleaseMemObject(pac_input_a_);
    clReleaseMemObject(pac_input_b_);
#endif
   // clReleaseProgram(program_);
   // clReleaseContext(context_);
  }
  cl_command_queue command_;
  cl_mem pac_input_a_;
  cl_mem pac_input_b_;
};

extern BWAOCLEnv* opencl_env;
#endif

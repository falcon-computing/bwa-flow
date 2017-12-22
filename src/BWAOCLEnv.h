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
    BWAOCLEnv(const char* bin_path, 
        const char* pac_path, 
        const char* kernel_name) :OpenCLEnv(bin_path, kernel_name) {
#ifdef INTEL_FPGA
      DLOG(ERROR) << "Intel OpenCL is not supported in this version";
#else
      FILE* pac_inp = fopen(pac_path, "rb");
      int pac_size = 0;
      fseek(pac_inp, 0, SEEK_END);
      pac_size = ftell(pac_inp);
      fseek(pac_inp, 0, SEEK_SET);
      pac_size = (pac_size+3)/4;

      int *pac = (int*)malloc(sizeof(int)*pac_size);
      fread(pac, sizeof(int), pac_size, pac_inp);
      fclose(pac_inp);

      // transfer PAC reference to all the devices
      for (int i = 0; i < device_envs_.size(); i++) {

        cl_context context = device_envs_[i].context;
        cl_command_queue cmd = device_envs_[i].cmd;

        int err = 0;
        cl_mem_ext_ptr_t ext_c, ext_d;
        ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
        ext_d.flags = XCL_MEM_DDR_BANK3; ext_d.obj = 0; ext_d.param = 0;

        pac_input_a_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
            sizeof(int)*pac_size, &ext_c, &err);
        OCL_CHECK(err, "fail to allocate ocl buffer for pac_a");

        pac_input_b_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
            sizeof(int)*pac_size, &ext_d, &err);
        OCL_CHECK(err, "fail to allocate ocl buffer for pac_b");

        cl_event event[2];
        OCL_CHECKRUN(clEnqueueWriteBuffer(cmd, pac_input_a_, CL_FALSE, 0, sizeof(int)*pac_size, pac, 0, NULL, &event[0]),
            "fail to write buffer for pac_a");
        OCL_CHECKRUN(clEnqueueWriteBuffer(cmd, pac_input_b_, CL_FALSE, 0, sizeof(int)*pac_size, pac, 0, NULL, &event[1]),
            "fail to write buffer for pac_b");
        OCL_CHECKRUN(clWaitForEvents(2, event), "fail to write buffers");

        clReleaseEvent(event[0]);
        clReleaseEvent(event[1]);
      }
      free(pac);
#endif
    }

    ~BWAOCLEnv(){
      if (pac_input_a_) clReleaseMemObject(pac_input_a_);
      if (pac_input_b_) clReleaseMemObject(pac_input_b_);
    }
    cl_mem pac_input_a_ = NULL;
    cl_mem pac_input_b_ = NULL;
};

extern BWAOCLEnv* opencl_env;
#endif

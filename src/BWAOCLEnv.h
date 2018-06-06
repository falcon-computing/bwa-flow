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
#include "config.h"
#include "OpenCLEnv.h"
#include "util.h"

// NOTE: the FPGA kernel requires a different layout for the pac
// buffer, hence the _set_pac() and _get_pac() macros are different
// from bwa/bntseq.c
// _get_pac_2() is introduced to get the origin pac array from
// aux->idx, reorganizing it for FPGA
#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((l&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((l&3)<<1)&3)
#define _get_pac_2(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

class BWAOCLEnv : public OpenCLEnv{
 public:
  BWAOCLEnv(
      const char* bin_path, 
      const char* kernel_name) :OpenCLEnv(bin_path, kernel_name) 
  {
    // get full pac array from pac
    char* pac;
    int64_t pac_size = get_full_pac(pac);

    int err = 0;
#ifdef XILINX_FPGA
    cl_mem_ext_ptr_t ext_c, ext_d;
    ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
    ext_d.flags = XCL_MEM_DDR_BANK3; ext_d.obj = 0; ext_d.param = 0;

    pac_input_a_ = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
        pac_size, &ext_c, &err);
    pac_input_b_ = clCreateBuffer(context_, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
        pac_size, &ext_d, &err);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("Failed to create reference OpenCL buffer!");
    }
    cl_event event[2];
    cl_command_queue command = getCmdQueue();
    err = clEnqueueWriteBuffer(command, pac_input_a_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[0]);
    err = clEnqueueWriteBuffer(command, pac_input_b_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[1]);
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
#ifdef XILINX_FPGA
    clReleaseMemObject(pac_input_a_);
    clReleaseMemObject(pac_input_b_);
#endif
  }

  static int64_t get_full_pac(char* &pac) { 
    int64_t l_pac = aux->idx->bns->l_pac;
    int64_t pac_size = (aux->idx->bns->l_pac * 2 + 3) / 4;

    pac = (char*)calloc(pac_size, 1);

    int64_t k = l_pac;
    
    // forward, using different set/get_pac schemes
    for (int64_t l = 0; l < l_pac; l++) {
      _set_pac(pac, l, _get_pac_2(aux->idx->pac, l));
    }
    // backward
    for (int64_t l = l_pac - 1; l >= 0; --l) {
      _set_pac(pac, k, 3-_get_pac(pac, l));
      k++;
    }

    return pac_size;
  }

  cl_mem pac_input_a_;
  cl_mem pac_input_b_;
};

extern BWAOCLEnv* opencl_env;
#endif

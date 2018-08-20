#ifndef BWAOCLENV_H
#define BWAOCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <glog/logging.h>
#include <string>
#include <vector>
#include <stdexcept>

#include <CL/opencl.h>

#include "kflow/Queue.h"
#include "config.h"
#include "OpenCLEnv.h"
#include "util.h"
#include "allocation_wrapper.h"

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
  BWAOCLEnv(const char* bin_path, 
            const char* kernel_name)
  : OpenCLEnv(bin_path, kernel_name) 
  {
    //initPAC();
    initBWT();
  }

  ~BWAOCLEnv() {
    //releasePAC();
    releaseBWT();
  }

  void initPAC() {
    // get full pac array from pac
    char* pac;
    int64_t pac_size = get_full_pac(pac);

#ifdef XILINX_FPGA
    cl_int err = 0;
    cl_mem_ext_ptr_t ext_c, ext_d;
    ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
    ext_d.flags = XCL_MEM_DDR_BANK3; ext_d.obj = 0; ext_d.param = 0;
    // transfer PAC reference to all the devices
    for (int i = 0; i < device_envs_.size(); i++) {
        cl_context context = device_envs_[i].context;
        cl_command_queue cmd = device_envs_[i].cmd;

        cl_mem pac_input_a_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
            pac_size, &ext_c, &err);
        cl_mem pac_input_b_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
            pac_size, &ext_d, &err);
        if (err != CL_SUCCESS) {
            throw std::runtime_error("Failed to create reference OpenCL buffer!");
        }
        cl_event event[2];
        err  = clEnqueueWriteBuffer(cmd, pac_input_a_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[0]);
        err |= clEnqueueWriteBuffer(cmd, pac_input_b_, CL_TRUE, 0, pac_size, pac, 0, NULL, &event[1]);
        clWaitForEvents(2, event);
        if (err != CL_SUCCESS) {
            throw std::runtime_error("Failed to write reference to DDR!");
        }
        clReleaseEvent(event[0]);
        clReleaseEvent(event[1]);
        pac_input_a_list_.push_back(pac_input_a_);
        pac_input_b_list_.push_back(pac_input_b_);
    }
#else
    DLOG(ERROR) << "This feature is currently only supported in Xilinx";
#endif
    free(pac);
  }

  void releasePAC() {
#ifdef XILINX_FPGA
    for (int i = 0; i < pac_input_a_list_.size(); i++) {
      clReleaseMemObject(pac_input_a_list_[i]);
    }
    for (int i = 0; i < pac_input_b_list_.size(); i++) {
      clReleaseMemObject(pac_input_b_list_[i]);
    }
#endif
  }

  void initBWT() {
    uint32_t *bwt          = aux->idx->bwt->bwt;
    size_t    bwt_size     = aux->idx->bwt->bwt_size*sizeof(uint32_t);
    uint64_t  bwt_param[6] = {aux->idx->bwt->primary,
                              aux->idx->bwt->L2[0],
                              aux->idx->bwt->L2[1],
                              aux->idx->bwt->L2[2],
                              aux->idx->bwt->L2[3],
                              aux->idx->bwt->L2[4]};

    cl_int err = 0;
    cl_mem_ext_ptr_t ext[4];
    ext[0].flags = XCL_MEM_DDR_BANK0; ext[0].obj = 0; ext[0].param=0;
    ext[1].flags = XCL_MEM_DDR_BANK1; ext[1].obj = 0; ext[1].param=0;
    ext[2].flags = XCL_MEM_DDR_BANK2; ext[2].obj = 0; ext[2].param=0;
    ext[3].flags = XCL_MEM_DDR_BANK3; ext[3].obj = 0; ext[3].param=0;

    // transfer BWT reference to all the devices
    for (int i = 0; i < device_envs_.size(); i++) {
      cl_context context = device_envs_[i].context;
      cl_command_queue cmd = device_envs_[i].cmd;

      cl_mem *bwt_ = new cl_mem[4];
      cl_mem *bwt_param_ = new cl_mem[4];
      for (int p = 0; p < 4; p++) {
        bwt_[p] = clCreateBuffer( context,
                                  CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                                  bwt_size,
                                  &(ext[ p ]),
                                  &err );
        if (err != CL_SUCCESS)
          throw std::runtime_error("Failed to create BWT reference OpenCL buffer!");
        bwt_param_[p] = clCreateBuffer( context,
                                        CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                                        6*sizeof(uint64_t),
                                        &(ext[ p ]),
                                        &err );
        if (err != CL_SUCCESS)
          throw std::runtime_error("Failed to create BWT reference (param) OpenCL buffer!");
      }

      cl_event bwt_events[4], bwt_param_events[4];
      for (int p = 0; p < 4; p++) {
        err |= clEnqueueWriteBuffer( cmd, bwt_[p], CL_FALSE, 0, bwt_size,
                                     bwt, 0, NULL, &bwt_events[p]);
        err |= clEnqueueWriteBuffer( cmd, bwt_param_[p], CL_FALSE, 0, 6*sizeof(uint64_t),
                                     bwt_param, 0, NULL, &bwt_param_events[p]);
      }
      err |= clWaitForEvents(4, bwt_events);
      err |= clWaitForEvents(4, bwt_param_events);
      if (err != CL_SUCCESS)
        throw std::runtime_error("Failed to write BWT reference to DDR");

      for (int p = 0; p < 4; p++) {
        clReleaseEvent(bwt_events[p]);
        clReleaseEvent(bwt_param_events[p]);
      }

      bwt_list_.push_back(bwt_);
      bwt_param_list_.push_back(bwt_param_);
    }
  }

  void releaseBWT() {
    for (int i = 0; i < bwt_list_.size(); i++)
      for (int p = 0; p < 4; p++)
        clReleaseMemObject(bwt_list_[i][p]);

    for (int i = 0; i < bwt_param_list_.size(); i++)
      for (int p = 0; p < 4; p++)
        clReleaseMemObject(bwt_param_list_[i][p]);
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

 public:
  std::vector<cl_mem*> bwt_list_;
  std::vector<cl_mem*> bwt_param_list_;
  std::vector<cl_mem> pac_input_a_list_;
  std::vector<cl_mem> pac_input_b_list_;
};

extern BWAOCLEnv* opencl_env;
#endif

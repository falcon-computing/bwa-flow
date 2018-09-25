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

struct cl_pe {
  int              pe_id;
  int              bank_id;
  std::string      type;
  cl_accx         *accx; 
  cl_command_queue cmd;
  cl_pe(): pe_id(-1), bank_id(-1), accx(NULL), cmd(NULL) {;}
};

const cl_pe NULL_PE;

class BWAOCLEnv : public OpenCLEnv{
 public:
  BWAOCLEnv()
    : OpenCLEnv(get_group_sizes(), get_bin_paths())
  {
    num_pe_ = 0;
    if (!FLAGS_no_use_sw_fpga)
      initPAC();
    if (!FLAGS_no_use_smem_fpga)
      initBWT();
  }

  ~BWAOCLEnv() {
    if (!FLAGS_no_use_sw_fpga)
      releasePAC();
    if (!FLAGS_no_use_smem_fpga)
      releaseBWT();
  }

  void initPAC() {
    sw_num_pe_ = 0;
    // get full pac array from pac
    char* pac;
    int64_t pac_size = get_full_pac(pac);

#ifdef XILINX_FPGA
    cl_int err = 0;
    cl_mem_ext_ptr_t ext_c, ext_d;
    ext_c.flags = XCL_MEM_DDR_BANK1; ext_c.obj = 0; ext_c.param = 0;
    ext_d.flags = XCL_MEM_DDR_BANK1; ext_d.obj = 0; ext_d.param = 0;
    // transfer PAC reference to all the devices
    for (int i = 0; i < device_envs_.size(); i++) {
      if (device_envs_[i].accx_group_id != 0) continue;

      cl_context context = device_envs_[i].context;
      cl_command_queue cmd = device_envs_[i].cmd;

      cl_mem pac_input_a_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
          pac_size, &ext_c, &err);
        if (err == CL_INVALID_CONTEXT) DLOG(ERROR) << "ctx";
        if (err == CL_INVALID_VALUE) DLOG(ERROR) << "val";
        if (err == CL_INVALID_BUFFER_SIZE) DLOG(ERROR) << "bsz";
        if (err == CL_INVALID_HOST_PTR) DLOG(ERROR) << "ptr";
        if (err == CL_MEM_OBJECT_ALLOCATION_FAILURE) DLOG(ERROR) << "alc";
        if (err == CL_OUT_OF_HOST_MEMORY) DLOG(ERROR) << "oom";
      cl_mem pac_input_b_ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
          pac_size, &ext_d, &err);
        if (err == CL_INVALID_CONTEXT) DLOG(ERROR) << "ctx";
        if (err == CL_INVALID_VALUE) DLOG(ERROR) << "val";
        if (err == CL_INVALID_BUFFER_SIZE) DLOG(ERROR) << "bsz";
        if (err == CL_INVALID_HOST_PTR) DLOG(ERROR) << "ptr";
        if (err == CL_MEM_OBJECT_ALLOCATION_FAILURE) DLOG(ERROR) << "alc";
        if (err == CL_OUT_OF_HOST_MEMORY) DLOG(ERROR) << "oom";
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

      cl_pe pe;
      pe.pe_id = sw_num_pe_++;
      pe.bank_id = 1;
      pe.type = "sw";
      pe.accx = &device_envs_[i];
      pe.cmd = clCreateCommandQueue(pe.accx->context, pe.accx->device_id, CL_QUEUE_PROFILING_ENABLE, &err);
      //pe.cmd = clCreateCommandQueue(pe.accx->context, pe.accx->device_id, 0, &err);
      OCL_CHECK(err, "failed to create cmd_queue");
      sw_pe_list_.push_back(pe);
    }
#else
    DLOG(ERROR) << "This feature is currently only supported in Xilinx";
#endif
    free(pac);
    sw_fpga_thread_ = sw_pe_list_.size();
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
    smem_num_pe_ = 0;
    uint32_t *bwt           = aux->idx->bwt->bwt;
    uint64_t  bwt_param[7]  = {aux->idx->bwt->primary,
                               aux->idx->bwt->L2[0],
                               aux->idx->bwt->L2[1],
                               aux->idx->bwt->L2[2],
                               aux->idx->bwt->L2[3],
                               aux->idx->bwt->L2[4],
                               (aux->idx->bwt->bwt_size+15)/16};
    size_t    bwt_bytes_pad = bwt_param[6]*16*sizeof(uint32_t);
    size_t    bwt_param_num = 7;

    cl_int err = 0;
    cl_mem_ext_ptr_t ext[4];
    ext[0].flags = XCL_MEM_DDR_BANK0; ext[0].obj = 0; ext[0].param=0;
    ext[1].flags = XCL_MEM_DDR_BANK1; ext[1].obj = 0; ext[1].param=0;
    ext[2].flags = XCL_MEM_DDR_BANK2; ext[2].obj = 0; ext[2].param=0;
    ext[3].flags = XCL_MEM_DDR_BANK3; ext[3].obj = 0; ext[3].param=0;

    // transfer BWT reference to all the devices
    for (int i = 0; i < device_envs_.size(); i++) {
      if (device_envs_[i].accx_group_id != 0) continue;

      cl_context context = device_envs_[i].context;
      cl_command_queue cmd = device_envs_[i].cmd;

      cl_mem bwt_input[4], bwt_param_input[4];
      for (int p = 0; p < 4; p++) {
        if (p == 1) continue;
        bwt_input[p] = clCreateBuffer( context,
                                  CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                                  bwt_bytes_pad, /*create the ocl buffer with padded size*/
                                  &(ext[ p ]),
                                  &err );
        if (err != CL_SUCCESS)
          throw std::runtime_error("Failed to create BWT reference OpenCL buffer!");
        bwt_param_input[p] = clCreateBuffer( context,
                                             CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX,
                                             bwt_param_num*sizeof(uint64_t),
                                             &(ext[ p ]),
                                             &err );
        if (err != CL_SUCCESS)
          throw std::runtime_error("Failed to create BWT reference (param) OpenCL buffer!");
      }

      cl_event bwt_events[4], bwt_param_events[4];
      for (int p = 0; p < 4; p++) {
        if (p == 1) continue;
        err |= clEnqueueWriteBuffer( cmd, bwt_input[p], CL_FALSE, 0, aux->idx->bwt->bwt_size*sizeof(uint32_t), /*transfer bwt with original size*/
                                     bwt, 0, NULL, &bwt_events[p]);
        err |= clEnqueueWriteBuffer( cmd, bwt_param_input[p], CL_FALSE, 0, bwt_param_num*sizeof(uint64_t),
                                     bwt_param, 0, NULL, &bwt_param_events[p]);
      }
      cl_event tmp1[3] = {bwt_events[0], bwt_events[2], bwt_events[3]};
      cl_event tmp2[3] = {bwt_param_events[0], bwt_param_events[2], bwt_param_events[3]};
      err |= clWaitForEvents(3, tmp1);
      err |= clWaitForEvents(3, tmp2);
      if (err != CL_SUCCESS)
        throw std::runtime_error("Failed to write BWT reference to DDR");

      for (int p = 0; p < 4; p++) {
        if (p == 1) continue;
        clReleaseEvent(bwt_events[p]);
        clReleaseEvent(bwt_param_events[p]);
      }

      for (int p = 0; p < 4; p++) {
        if (p == 1) continue;
        bwt_list_.push_back(bwt_input[p]);
        bwt_param_list_.push_back(bwt_param_input[p]);
      }

      for (int bank_id = 0; bank_id < 4; bank_id++) {
        if (bank_id == 1) continue;
        cl_pe pe;
        pe.pe_id = smem_num_pe_++;
        pe.bank_id = bank_id;
        pe.type = "smem";
        pe.accx = &device_envs_[i];
        pe.cmd = clCreateCommandQueue(pe.accx->context, pe.accx->device_id, CL_QUEUE_PROFILING_ENABLE|CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
        //pe.cmd = clCreateCommandQueue(pe.accx->context, pe.accx->device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
        OCL_CHECK(err, "failed to create cmd_queue");
        smem_pe_list_.push_back(pe);
      }
    }
    smem_fpga_thread_ = smem_pe_list_.size();
  }

  void releaseBWT() {
    for (int i = 0; i < bwt_list_.size(); i++)
        clReleaseMemObject(bwt_list_[i]);

    for (int i = 0; i < bwt_param_list_.size(); i++)
        clReleaseMemObject(bwt_param_list_[i]);
  }

  cl_pe getPE(std::string type) {
    boost::lock_guard<OpenCLEnv> guard(*this);
    uint32_t tid = getTid();
    DLOG(INFO) << tid;
    if (pe_registry_.count(tid)) {
      pe_refs_[tid]++;
      return pe_registry_[tid];
    }
    else {
      std::vector<cl_pe> *pe_list = NULL;
      DLOG(INFO) << "type: " <<  type;
      if (type == "sw") pe_list = &sw_pe_list_;
      else if (type == "smem") pe_list = &smem_pe_list_;

      if (pe_list && !pe_list->empty()) {
        cl_pe pe = pe_list->back();
        pe_list->pop_back();
        pe_registry_[tid] = pe;
        pe_refs_[tid] = 0;
        return pe;
      }
      else {
        DLOG(ERROR) << "No more " << type << " pe available on this platform";
        return NULL_PE;
      }
    }
  }

  void releasePE(cl_pe pe) {
    boost::lock_guard<OpenCLEnv> guard(*this);
    uint32_t tid = getTid();
    if (--pe_refs_[tid] == 0) {
      pe_refs_.erase(tid);
      pe_registry_.erase(tid);
      if (pe.type == "sw") sw_pe_list_.push_back(pe);
      else if (pe.type == "smem") smem_pe_list_.push_back(pe);
    }
  }

 private:
  std::vector<int> &get_group_sizes() {
    group_sizes_ = new std::vector<int>();
    
    int max_devices = getMaxNumDevices();
    if (FLAGS_max_fpga_thread >= 0) {
      max_devices = std::min(FLAGS_max_fpga_thread, int(max_devices));
    }

    group_sizes_->push_back(max_devices);

    return *group_sizes_;
  }

  std::vector<const char*> &get_bin_paths() {
    bin_paths_ = new std::vector<const char*>();
    bin_paths_->push_back(FLAGS_fpga_path.c_str());

    return *bin_paths_;
  }

 public:
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

  // get current thread id
  // using the same code from googlelog/src/utilities.cc
  // without OS checking
 private:
  uint32_t getTid() {
    static bool lacks_gettid = false;

    if (!lacks_gettid) {
      pid_t tid = syscall(__NR_gettid);
      if (tid != -1) {
        return (uint32_t)tid;
      }
      // Technically, this variable has to be volatile, but there is a small
      // performance penalty in accessing volatile variables and there should
      // not be any serious adverse effect if a thread does not immediately see
      // the value change to "true".
      lacks_gettid = true;
    }

    // If gettid() could not be used, we use one of the following.
    return (uint32_t)getpid();
  }

 public:
  std::vector<cl_mem> bwt_list_;
  std::vector<cl_mem> bwt_param_list_;
  std::vector<cl_mem> pac_input_a_list_;
  std::vector<cl_mem> pac_input_b_list_;

  std::vector<int> *group_sizes_;
  std::vector<const char*> *bin_paths_;

  int num_pe_;
  int sw_num_pe_;
  int smem_num_pe_;
  std::vector<cl_pe> sw_pe_list_;
  std::vector<cl_pe> smem_pe_list_;

  int sw_fpga_thread_;
  int smem_fpga_thread_;

  std::map<uint32_t, cl_pe> pe_registry_;
  std::map<uint32_t, int>   pe_refs_;
};

extern BWAOCLEnv* opencl_env;
#endif

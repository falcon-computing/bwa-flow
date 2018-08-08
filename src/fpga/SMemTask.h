#ifndef SMEM_TASK_H
#define SMEM_TASK_H

#include "bwa_wrapper.h"
#include "Task.h"
#include "FPGAAgent.h"
#include "BWAOCLEnv.h"

class SMemTask : public Task {

 public:
  SMemTask(BWAOCLEnv* env);
  ~SMemTask();

  void start(SMemTask* prev_task);
  void finish();

  cl_mem     i_seq_buf;
  cl_mem     o_mem_buf;
  cl_mem     o_num_buf;

  size_t     i_seq_size;
  size_t     o_mem_size;
  size_t     o_num_size;

  char      *i_seq_data;
  bwtintv_t *o_mem_data;
  int       *o_num_data;

  int        i_seq_base_idx;
  int        i_seq_num;

  size_t     max_i_seq_num_ = 320;
  size_t     max_i_seq_len_ = 150;
  size_t     max_intv_alloc_ = 200;

 private:
  FPGAAgent* agent_;
};
#endif

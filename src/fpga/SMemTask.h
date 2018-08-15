#ifndef SMEM_TASK_H
#define SMEM_TASK_H

#include "bwa_wrapper.h"
#include "Task.h"
#include "FPGAAgent.h"
#include "BWAOCLEnv.h"

#define SMEM_BANK_NUM 4

class SMemTask : public Task {

 public:
  SMemTask(BWAOCLEnv* env);
  ~SMemTask();

  void start(SMemTask* prev_task);
  void finish();

  cl_mem     i_seq_buf[SMEM_BANK_NUM];
  cl_mem     i_seq_len_buf[SMEM_BANK_NUM];
  cl_mem     o_mem_buf[SMEM_BANK_NUM];
  cl_mem     o_num_buf[SMEM_BANK_NUM];

  size_t     i_seq_size[SMEM_BANK_NUM];
  size_t     i_seq_len_size[SMEM_BANK_NUM];
  size_t     o_mem_size[SMEM_BANK_NUM];
  size_t     o_num_size[SMEM_BANK_NUM];

  uint8_t   *i_seq_data[SMEM_BANK_NUM];
  uint8_t   *i_seq_len_data[SMEM_BANK_NUM];
  bwtintv_t *o_mem_data[SMEM_BANK_NUM];
  int       *o_num_data[SMEM_BANK_NUM];

  int        i_seq_base_idx;
  int        i_seq_num[SMEM_BANK_NUM];

  size_t     max_i_seq_num_ = 1024;
  size_t     max_i_seq_len_ = 256;
  size_t     max_intv_alloc_ = 256;

 private:
  FPGAAgent* agent_;
};
#endif

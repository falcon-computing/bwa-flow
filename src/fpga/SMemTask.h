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
  void start(SMemTask* prev_task, uint64_t &write_ts, uint64_t &enq_ts);
  void finish();
  void finish(uint64_t &deq_ts, uint64_t &read_ts);

  cl_mem     i_bwt;
  cl_mem     i_bwt_param;

  cl_mem     i_seq_buf;
  cl_mem     i_seq_len_buf;
  cl_mem     o_mem_buf;
  cl_mem     o_num_buf;

  size_t     i_seq_size;
  size_t     i_seq_len_size;
  size_t     o_mem_size;
  size_t     o_num_size;

  uint8_t   *i_seq_data;
  uint8_t   *i_seq_len_data;
  bwtintv_t *o_mem_data;
  int       *o_num_data;

  int        i_seq_base_idx;
  int        i_seq_num;

  static const size_t max_i_seq_num_;
  static const size_t max_i_seq_len_;
  static const size_t max_intv_alloc_;

 private:
  FPGAAgent* agent_;
};
#endif

#ifndef SW_TASK_H
#define SW_TASK_H

#include "bwa_wrapper.h"
#include "FPGAAgent.h"
#include "BWAOCLEnv.h"

class SWTask {

 public:
  SWTask(BWAOCLEnv* env, int chunk_size);
  ~SWTask();

  void start(SWTask* prev_task);
  void finish();

  int     i_size[2];
  int     o_size[2];
  cl_mem  i_buf[2];
  cl_mem  o_buf[2];
  char*   i_data[2];
  short*  o_data[2];

  cl_event*       events;
  mem_alnreg_t**  region_batch;
  mem_chain_t**   chain_batch;

  size_t     max_i_size_;
  size_t     max_o_size_;
 private:
  FPGAAgent* agent_;
};
#endif

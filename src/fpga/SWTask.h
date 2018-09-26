#ifndef SW_TASK_H
#define SW_TASK_H

#include "bwa_wrapper.h"
#include "Task.h"
#include "FPGAAgent.h"
#include "BWAOCLEnv.h"

class SWTask : public Task {

 public:
  SWTask(BWAOCLEnv* env, int chunk_size);
  ~SWTask();

  void start(SWTask* prev_task);
  void finish();
  void redo();

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

  int        start_seq = 0;
  int        end_seq = 0;

 private:
  FPGAAgent* agent_;

  boost::thread helper_;
  boost::mutex mtx_;
  boost::atomic<int> state_;
  boost::atomic<SWTask*> prv_task_;

 private:
  void start_func(SWTask* prev_task);
  void finish_func();
  void redo_func();

  void helper_func();

};
#endif

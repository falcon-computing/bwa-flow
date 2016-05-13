#ifndef FPGAAGENT_H
#define FPGAAGENT_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>

#include "OpenCLEnv.h"

#define FPGA_RET_PARAM_NUM 5

class FPGAAgent {
 public:
  FPGAAgent(OpenCLEnv* env, 
      int chunk_size,
      uint64_t buf_size = 32*1024*1024);

  ~FPGAAgent();

  void writeInput(void* host_ptr, uint64_t size, int cnt);
  void readOutput(void* host_ptr, uint64_t size, int cnt);
  void start(int task_num, int cnt);

 private:
  OpenCLEnv*     env_;
  const uint64_t max_buf_size_;
  const int      chunk_size_;
  cl_mem         input_buf_[2];
  cl_mem         output_buf_[2];
  void*          host_buf_[2];
};

#endif

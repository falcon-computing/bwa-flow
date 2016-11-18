#ifndef FPGAAGENT_H
#define FPGAAGENT_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>

#include "OpenCLEnv.h"

#define FPGA_RET_PARAM_NUM 5

extern OpenCLEnv* opencl_env;

class FPGAAgent {
 public:
  FPGAAgent(OpenCLEnv* env, 
      int chunk_size,
      uint64_t buf_size = 32*1024*1024);

  ~FPGAAgent();

  void writeInput(void* host_ptr, uint64_t size, int cnt, int bank);
  void readOutput(void* host_ptr, uint64_t size, int cnt, int bank);
  void start(int size_a, int size_b, int cnt);
  void wait(int cnt);
  bool pending(int cnt);

 private:
  OpenCLEnv*     env_;
  const uint64_t max_buf_size_;
  const int      chunk_size_;
  cl_mem         input_buf_a[2];
  cl_mem         input_buf_b[2];
  cl_mem         output_buf_a[2];
  cl_mem         output_buf_b[2];
  FPGATask*      fpga_task_[2];
};

#endif

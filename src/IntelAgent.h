#ifndef INTELAGENT_H
#define INTELAGENT_H
#include "FPGAAgent.h"
#include "SWTask.h"

#define AOCL_ALIGNMENT 64
inline void *sw_malloc(size_t size, int data_width) {
  size_t aligned_size = data_width*((size/AOCL_ALIGNMENT)+1)*AOCL_ALIGNMENT;
  void *result = NULL;
  posix_memalign(&result, AOCL_ALIGNMENT, aligned_size);
  return result;
}

class IntelAgent : public FPGAAgent {
 public:
  IntelAgent(OpenCLEnv* env, SWTask* task);
  ~IntelAgent();

  void writeInput(void* host_ptr, int size, int bank);
  void readOutput(void* host_ptr, int size, int bank);
  void start(FPGAAgent* prev_agent = NULL);
  void finish();

 private:
  SWTask*          task_;
  cl_kernel        kernels_[4];
  cl_command_queue cmd_[4];
  cl_mem           i_buf_[2];
  cl_mem           o_buf_[2];
  cl_event         events_[4];
};

#endif

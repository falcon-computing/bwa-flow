#ifndef XCL_AGENT_H
#define XCL_AGENT_H
#include "FPGAAgent.h"
#include "SWTask.h"

#define XCL_ALIGNMENT 4096
inline void *sw_malloc(size_t size, int data_width) {
  size_t aligned_size = data_width*((size/XCL_ALIGNMENT)+1)*XCL_ALIGNMENT;
  void *result = NULL;
  posix_memalign(&result, XCL_ALIGNMENT, aligned_size);
  return result;
}

class XCLAgent : public FPGAAgent {
 public:
  XCLAgent(OpenCLEnv* env, SWTask* task);
  ~XCLAgent();

  void writeInput(void* host_ptr, int size, int bank);
  void readOutput(void* host_ptr, int size, int bank);
  void start(FPGAAgent* prev_agent = NULL);
  void finish();

 private:
  cl_context       context_;
  cl_kernel        kernel_;
  cl_command_queue cmd_;
  cl_mem           i_buf_[2];
  cl_mem           o_buf_[2];
  
  cl_event         kernel_event_;
  cl_event         write_events_[2];

  SWTask*          task_;
};

#endif

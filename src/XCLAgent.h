#ifndef XCL_AGENT_H
#define XCL_AGENT_H
#include "FPGAAgent.h"

class XCLAgent : public FPGAAgent {
 public:
  XCLAgent(BWAOCLEnv* env, SWTask* task);
  ~XCLAgent();

  void writeInput(cl_mem buf, void* host_ptr, int size, int bank);
  void readOutput(cl_mem buf, void* host_ptr, int size, int bank);
  void start(SWTask* task, FPGAAgent* prev_agent = NULL);
  void finish();
  void fence();

 private:
  cl_kernel        kernel_;
  cl_command_queue cmd_;
  cl_event         kernel_event_;
  cl_event         write_events_[2];
  BWAOCLEnv*       env_;
  bool             valid_2nd_event_;
};

#endif

#ifndef XCL_AGENT_H
#define XCL_AGENT_H

#include "FPGAAgent.h"
#include "SWTask.h"

class XCLAgent : public FPGAAgent {
 public:
  XCLAgent(BWAOCLEnv* env, SWTask* task);
  ~XCLAgent();

  void writeInput(cl_mem buf, void* host_ptr, int size, int bank);
  void readOutput(cl_mem buf, void* host_ptr, int size, int bank);
  void start(Task* task, FPGAAgent* prev_agent = NULL);
  void finish();
  void fence();

 private:
  cl_pe         pe_;
  BWAOCLEnv*    env_;

  cl_kernel     kernel_;

  cl_event      kernel_event_;
  cl_event      write_events_;
  cl_event      read_events_;

  uint64_t      kernel_time_;
  uint64_t      kernel_invks_;
  uint64_t      writing_time_;
  uint64_t      reading_time_;
};

#endif

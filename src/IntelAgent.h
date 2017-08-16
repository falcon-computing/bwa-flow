#ifndef INTELAGENT_H
#define INTELAGENT_H
#include "FPGAAgent.h"

class IntelAgent : public FPGAAgent {
 public:
  IntelAgent(OpenCLEnv* env);
  ~IntelAgent();

  void writeInput(cl_mem buf, void* host_ptr, int size, int bank);
  void readOutput(cl_mem buf, void* host_ptr, int size, int bank);
  void start(SWTask* task, FPGAAgent* prev_agent = NULL);
  void finish();

 private:
  cl_kernel        kernels_[4];
  cl_command_queue cmd_[4];
  cl_event         events_[4];
};

#endif

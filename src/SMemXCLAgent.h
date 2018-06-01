#ifndef SMEM_XCL_AGENT_H
#define SMEM_XCL_AGENT_H
#include "FPGAAgent.h"
#include "SMemTask.h"
#include <vector>

class SMemXCLAgent : public FPGAAgent {
 public:
  SMemXCLAgent(BWAOCLEnv* env, SMemTask* task);
  ~SMemXCLAgent();

  void writeInput(cl_mem buf, void* host_ptr, int size, int bank);
  void readOutput(cl_mem buf, void* host_ptr, int size, int bank);
  void start(Task* task, FPGAAgent* prev_agent = NULL);
  void finish();

  void createBuffer(SMemTask* task);
  void releaseBuffer(SMemTask* task);

 private:
  cl_kernel        kernel_;
  cl_event         kernel_event_;
  cl_event         write_events_;
  BWAOCLEnv*       env_;
};

#endif

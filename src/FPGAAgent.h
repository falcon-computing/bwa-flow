#ifndef FPGAAGENT_H
#define FPGAAGENT_H
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>

#include "bwa_wrapper.h"
#include "BWAOCLEnv.h"

class SWTask;

class FPGAAgent {
 public:
  virtual void writeInput(cl_mem buf, void* host_ptr, int size, int bank) = 0;
  virtual void readOutput(cl_mem buf, void* host_ptr, int size, int bank) = 0;
  virtual void start(SWTask* task, FPGAAgent* prev_agent = NULL) = 0;
  virtual void finish() = 0;
};

#endif

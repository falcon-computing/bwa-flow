#ifndef FPGAAGENT_H
#define FPGAAGENT_H
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>

#include "bwa_wrapper.h"
#include "OpenCLEnv.h"

class SWTask;

class FPGAAgent {
 public:
  virtual void writeInput(void* host_ptr, int size, int bank) = 0;
  virtual void readOutput(void* host_ptr, int size, int bank) = 0;
  virtual void start(FPGAAgent* prev_agent = NULL) = 0;
  virtual void finish() = 0;
};

#endif

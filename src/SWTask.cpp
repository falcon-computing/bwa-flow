#include <unistd.h>

#include "bwa_wrapper.h"
#include "config.h"
#include "OpenCLEnv.h"
#include "SWTask.h"

#ifdef INTEL_FPGA
#include "IntelAgent.h"
#elif XILINX_FPGA
#include "XCLAgent.h"
#endif

SWTask::SWTask(OpenCLEnv* env, int chunk_size) {
  max_i_size_ = 32*1024*1024;
  max_o_size_ = 2*chunk_size*FPGA_RET_PARAM_NUM;

#ifdef INTEL_FPGA
  agent_ = new IntelAgent(env, this);
#elif XILINX_FPGA
  agent_ = new XCLAgent(env, this);
#endif

  for (int k = 0; k < 2; k++) {
    i_size[k] = 0;
    o_size[k] = 0;

    i_data[k] = (char*) sw_malloc(max_i_size_, sizeof(char));
    o_data[k] = (short*)sw_malloc(max_o_size_, sizeof(short));
  }
  region_batch = new mem_alnreg_t*[2*chunk_size];
  chain_batch  = new mem_chain_t*[2*chunk_size];
}

SWTask::~SWTask() {
  delete agent_;
}

void SWTask::start(SWTask* prev_task) {
  if (i_size[0]*sizeof(int) >= max_i_size_ || 
      i_size[1]*sizeof(int) >= max_i_size_ ||
      o_size[0] + o_size[1] >= max_o_size_) 
  {
    DLOG(ERROR) << "exceeding max memory size";
    throw std::runtime_error("exceeding max memory size");
  }
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Task info: " 
    << "i_size[0] = " << i_size[0] << ", "
    << "i_size[1] = " << i_size[1];
                               

  // write input
  uint64_t start_ts = getUs();
  agent_->writeInput(i_data[0], i_size[0]*sizeof(int), 0);
  agent_->writeInput(i_data[1], i_size[1]*sizeof(int), 1);

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Prepare input takes " << 
                               getUs() - start_ts << " us";

  if (prev_task->i_size[0] == 0 && prev_task->i_size[1] == 0) {
    agent_->start(NULL);
  }
  else {
    agent_->start(prev_task->agent_);
  }
}

void SWTask::finish() {
  uint64_t start_ts = getUs();

  agent_->readOutput(o_data[0], o_size[0]*sizeof(int), 0);
  agent_->readOutput(o_data[1], o_size[1]*sizeof(int), 1);

  // release events
  agent_->finish();

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Read OpenCL buffer takes " 
                               << getUs() - start_ts << " us";
}

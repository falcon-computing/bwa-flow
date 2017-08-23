#include <unistd.h>

#include "bwa_wrapper.h"
#include "config.h"
#include "BWAOCLEnv.h"
#include "SWTask.h"

#ifdef INTEL_FPGA
#include "IntelAgent.h"
#elif XILINX_FPGA
#include "XCLAgent.h"
#endif

#ifdef INTEL_FPGA
#define AOCL_ALIGNMENT 64
inline void *sw_malloc(size_t size, int data_width) {
  size_t aligned_size = data_width*((size/AOCL_ALIGNMENT)+1)*AOCL_ALIGNMENT;
  void *result = NULL;
  posix_memalign(&result, AOCL_ALIGNMENT, aligned_size);
  return result;
}
#else
inline void *sw_malloc(size_t size, int data_width) {
  return malloc(size*data_width);
}
#endif


SWTask::SWTask(BWAOCLEnv* env, int chunk_size) {
  max_i_size_ = 32*1024*1024;
  max_o_size_ = 2*chunk_size*FPGA_RET_PARAM_NUM;
#ifdef INTEL_FPGA
  agent_ = new IntelAgent(env, this);
#elif XILINX_FPGA
  agent_ = new XCLAgent(env, this);
#endif

//#ifdef XILINX_FPGA
  for (int k = 0; k < 2; k++) {
    i_size[k] = 0;
    o_size[k] = 0;

    i_data[k] = (char*) sw_malloc(max_i_size_, sizeof(char));
    o_data[k] = (short*)sw_malloc(max_o_size_, sizeof(short));
  }
//#elif INTEL_FPGA
//  for (int k = 0; k < 2; k++) {
//    i_size[k] = 0;
//    o_size[k] = 0;
//
//    i_data[k] = (char*) sw_malloc(max_i_size_, sizeof(char));
//    o_data[k] = (short*)sw_malloc(max_o_size_, sizeof(short));
//  }
//#endif
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
                               

  uint64_t start_ts = getUs();
  agent_->writeInput(i_buf[0], i_data[0], i_size[0]*sizeof(int), 0);
  agent_->writeInput(i_buf[1], i_data[1], i_size[1]*sizeof(int), 1);
  if (prev_task->i_size[0] == 0 && prev_task->i_size[1] == 0) {
    agent_->start(this, NULL);
  }
  else {
    agent_->start(this, prev_task->agent_);
  }
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Write OpenCL buffer takes " << 
                               getUs() - start_ts << " us";
}

void SWTask::finish() {
  uint64_t start_ts = getUs();

  agent_->readOutput(o_buf[0], o_data[0], o_size[0]*sizeof(int), 0);
  agent_->readOutput(o_buf[1], o_data[1], o_size[1]*sizeof(int), 1);

  // release events
  agent_->finish();

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Read OpenCL buffer takes " 
                               << getUs() - start_ts << " us";
}

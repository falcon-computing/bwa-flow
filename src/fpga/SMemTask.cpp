#include <unistd.h>

#include "bwa_wrapper.h"
#include "config.h"
#include "BWAOCLEnv.h"
#include "SMemTask.h"

#ifdef XILINX_FPGA
#include "SMemXCLAgent.h"
#endif

#ifdef INTEL_FPGA
#define AOCL_ALIGNMENT 64
inline void *smem_malloc(size_t size, int data_width) {
  size_t aligned_size = data_width*((size/AOCL_ALIGNMENT)+1)*AOCL_ALIGNMENT;
  void *result = NULL;
  posix_memalign(&result, AOCL_ALIGNMENT, aligned_size);
  return result;
}
#else
inline void *smem_malloc(size_t size, int data_width) {
  return malloc(size*data_width);
}
#endif

const size_t SMemTask::max_i_seq_num_ = 1024;
const size_t SMemTask::max_i_seq_len_ = 256;
const size_t SMemTask::max_intv_alloc_ = 256;

SMemTask::SMemTask(BWAOCLEnv* env) {
  i_seq_base_idx = -1;

  i_seq_num = 0;

  i_seq_data = (uint8_t *)smem_malloc(max_i_seq_len_*max_i_seq_num_, sizeof(uint8_t));
  i_seq_size = max_i_seq_len_*max_i_seq_num_ * sizeof(uint8_t);
  i_seq_len_data = (uint8_t *)smem_malloc(max_i_seq_num_, sizeof(uint8_t));
  i_seq_len_size = max_i_seq_num_ * sizeof(uint8_t);

  o_mem_data = (bwtintv_t *)smem_malloc(max_intv_alloc_*max_i_seq_num_, sizeof(bwtintv_t));
  o_mem_size = max_intv_alloc_*max_i_seq_num_ * sizeof(bwtintv_t);
  o_num_data = (int *)smem_malloc(max_i_seq_num_, sizeof(int));
  o_num_size = max_i_seq_num_ * sizeof(int);
  
#ifdef XILINX_FPGA
  agent_ = new SMemXCLAgent(env, this);

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Creating Buffer";
  ((SMemXCLAgent*)agent_)->createBuffer(this);
#endif

}

SMemTask::~SMemTask() {

#ifdef XILINX_FPGA
  ((SMemXCLAgent *)agent_)->releaseBuffer(this);
  delete agent_;
#endif

  free(i_seq_data);
  free(i_seq_len_data);
  free(o_mem_data);
  free(o_num_data);

}

void SMemTask::start(SMemTask* prev_task) {

#ifdef XILINX_FPGA
  uint64_t start_ts = getUs();
  //DLOG_IF(INFO, VLOG_IS_ON(3)) << "Try to write bank " << i;
  agent_->writeInput(i_seq_buf, i_seq_data, i_seq_size, 0);
  agent_->writeInput(i_seq_len_buf, i_seq_len_data, i_seq_len_size, 1);
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Writing OpenCL buffer takes "
                               << getUs() - start_ts << " us";
  
  start_ts = getUs();
  if (prev_task->i_seq_num == 0)
    agent_->start(this, NULL);
  else
    agent_->start(this, prev_task->agent_);
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Enqueuing tasks takes "
                               << getUs() - start_ts << " us";
#endif
}

void SMemTask::start(SMemTask* prev_task, uint64_t &write_ts, uint64_t &enq_ts) {
#ifdef XILINX_FPGA

  uint64_t start_ts = getUs();
  //DLOG_IF(INFO, VLOG_IS_ON(3)) << "Try to write bank " << i;
  agent_->writeInput(i_seq_buf, i_seq_data, i_seq_size, 0);
  agent_->writeInput(i_seq_len_buf, i_seq_len_data, i_seq_len_size, 1);
  write_ts += getUs() - start_ts;
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Writing OpenCL buffer takes "
                               << getUs() - start_ts << " us";
  
  start_ts = getUs();
  if (prev_task->i_seq_num == 0)
    agent_->start(this, NULL);
  else
    agent_->start(this, prev_task->agent_);
  enq_ts += getUs() - start_ts;
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Enqueuing tasks takes "
                               << getUs() - start_ts << " us";
#endif
}

void SMemTask::finish() {
#ifdef XILINX_FPGA
  uint64_t start_ts = getUs();
  agent_->readOutput(o_mem_buf, o_mem_data, o_mem_size, 0);
  agent_->readOutput(o_num_buf, o_num_data, o_num_size, 0);
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Reading OpenCL buffer takes " 
                               << getUs() - start_ts << " us";

  // release events
  agent_->finish();
#endif
}

void SMemTask::finish(uint64_t &deq_ts, uint64_t &read_ts) {
#ifdef XILINX_FPGA
  uint64_t start_ts = getUs();
  ((SMemXCLAgent*)agent_)->wait();
  deq_ts += getUs() - start_ts;
  start_ts = getUs();
  agent_->readOutput(o_mem_buf, o_mem_data, o_mem_size, 0);
  agent_->readOutput(o_num_buf, o_num_data, o_num_size, 0);
  read_ts += getUs() - start_ts;
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Reading OpenCL buffer takes " 
                               << getUs() - start_ts << " us";

  // release events
  agent_->finish();
#endif
}

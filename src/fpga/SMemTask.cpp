#include <unistd.h>

#include "bwa_wrapper.h"
#include "config.h"
#include "BWAOCLEnv.h"
#include "SMemTask.h"

#ifdef INTEL_FPGA
#include "SMemIntelAgent.h"
#elif XILINX_FPGA
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

SMemTask::SMemTask(BWAOCLEnv* env) {
  i_seq_base_idx = -1;

  for (int i=0; i<SMEM_BANK_NUM; i++) {
    i_seq_num[i] = 0;

    i_seq_data[i] = (uint8_t *)smem_malloc(max_i_seq_len_*max_i_seq_num_, sizeof(uint8_t));
    i_seq_size[i] = max_i_seq_len_*max_i_seq_num_ * sizeof(uint8_t);
    i_seq_len_data[i] = (uint8_t *)smem_malloc(max_i_seq_num_, sizeof(uint8_t));
    i_seq_len_size[i] = max_i_seq_num_ * sizeof(uint8_t);

    o_mem_data[i] = (bwtintv_t *)smem_malloc(max_intv_alloc_*max_i_seq_num_, sizeof(bwtintv_t));
    o_mem_size[i] = max_intv_alloc_*max_i_seq_num_ * sizeof(bwtintv_t);
    o_num_data[i] = (int *)smem_malloc(max_i_seq_num_, sizeof(int));
    o_num_size[i] = max_i_seq_num_ * sizeof(int);
  }
  total_i_seq_num = 0;
  
#ifdef INTEL_FPGA
  agent_ = new SMemIntelAgent(env, this);
#elif XILINX_FPGA
  agent_ = new SMemXCLAgent(env, this);
#endif

  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Creating Buffer";
  ((SMemXCLAgent*)agent_)->createBuffer(this);

}

SMemTask::~SMemTask() {

  ((SMemXCLAgent *)agent_)->releaseBuffer(this);

  for (int i=0; i<SMEM_BANK_NUM; i++) {
    free(i_seq_data[i]);
    free(i_seq_len_data[i]);
    free(o_mem_data[i]);
    free(o_num_data[i]);
  }

  delete agent_;
}

void SMemTask::start(SMemTask* prev_task) {

  uint64_t start_ts = getUs();
  for (int i=0; i<SMEM_BANK_NUM; i++) {
    //DLOG_IF(INFO, VLOG_IS_ON(3)) << "Try to write bank " << i;
    agent_->writeInput(i_seq_buf[i], i_seq_data[i], i_seq_size[i], i);
    agent_->writeInput(i_seq_len_buf[i], i_seq_len_data[i], i_seq_len_size[i], i);
  }
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Writing OpenCL buffer takes "
                               << getUs() - start_ts << " us";
  
  start_ts = getUs();
  
  if (prev_task->total_i_seq_num == 0)
    agent_->start(this, NULL);
  else
    agent_->start(this, prev_task->agent_);
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Enqueuing tasks takes "
                               << getUs() - start_ts << " us";
}

void SMemTask::finish() {
  uint64_t start_ts = getUs();
  for (int i=0; i<SMEM_BANK_NUM; i++) {
    agent_->readOutput(o_mem_buf[i], o_mem_data[i], o_mem_size[i], i);
    agent_->readOutput(o_num_buf[i], o_num_data[i], o_num_size[i], i);
  }
  DLOG_IF(INFO, VLOG_IS_ON(3)) << "Reading OpenCL buffer takes " 
                               << getUs() - start_ts << " us";

  // release events
  agent_->finish();
}

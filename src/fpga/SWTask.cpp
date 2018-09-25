#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#include "bwa_wrapper.h"
#include "config.h"
#include "allocation_wrapper.h"
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
  if (NULL == result) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  return result;
}
#else
inline void *sw_malloc(size_t size, int data_width) {
  return malloc(size*data_width);
}
#endif


SWTask::SWTask(BWAOCLEnv* env, int chunk_size) {
  max_i_size_ = 32*1024*1024;
  max_o_size_ = 2*2*chunk_size*FPGA_RET_PARAM_NUM;
  for (int k = 0; k < 2; k++) {
    i_size[k] = 0;
    o_size[k] = 0;

    i_data[k] = (char*) sw_malloc(max_i_size_, sizeof(char));
    o_data[k] = (short*)sw_malloc(max_o_size_, sizeof(short));
  }
  region_batch = new mem_alnreg_t*[2*chunk_size];
  if (NULL == region_batch) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  chain_batch  = new mem_chain_t*[2*chunk_size];
  if (NULL == chain_batch) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
#ifdef INTEL_FPGA
  agent_ = new IntelAgent(env, this);
#elif XILINX_FPGA
  agent_ = new XCLAgent(env, this);
#endif
  if (NULL == agent_) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }

  mtx_.lock();
  state_.store(0);
  helper_ = boost::thread(boost::bind(&SWTask::helper_func, this));
}

SWTask::~SWTask() {
  helper_.interrupt();
  //pthread_kill(helper_.native_handle(), 9);

  delete region_batch;
  delete chain_batch;
  for (int k = 0; k < 2; k++) {
    clReleaseMemObject(i_buf[k]);
    clReleaseMemObject(o_buf[k]);
    free(i_data[k]);
    free(o_data[k]);
  }
  delete agent_;
}

void SWTask::start(SWTask* prev_task) {
  state_.store(1);
  prv_task_.store(prev_task);
  mtx_.unlock();

  int wait = 0;
  while (!mtx_.try_lock()) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (wait++ >= 12000000)
      throw fpgaHangError("smithwater kernel stuck at start on fpga");
  }
}

void SWTask::start_func(SWTask* prev_task) {
  if (i_size[0]*sizeof(int) >= max_i_size_ || 
      i_size[1]*sizeof(int) >= max_i_size_ ||
      o_size[0] + o_size[1] >= max_o_size_) 
  {
    DLOG(ERROR) << "exceeding max memory size";
    throw std::runtime_error("exceeding max memory size");
  }
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Task info: " 
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

  agent_->readOutput(o_buf[0], o_data[0], o_size[0]*sizeof(int), 0);
  agent_->readOutput(o_buf[1], o_data[1], o_size[1]*sizeof(int), 1);

  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Enqueue write, kernel & read takes " <<
    getUs() - start_ts << " us";
}

void SWTask::finish() {
  state_.store(2);
  mtx_.unlock();

  int wait = 0;
  while (!mtx_.try_lock()) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (wait++ >= 12000000)
      throw fpgaHangError("smithwater kernel stuck at finish on fpga");
  }
}

void SWTask::finish_func() {
  uint64_t start_ts = getUs();
  agent_->finish();
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Dequeue write, kernel & read takes "
                               << getUs() - start_ts << " us";
}

void SWTask::redo() {
  state_.store(3);
  mtx_.unlock();

  int wait = 0;
  while (!mtx_.try_lock()) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (wait++ >= 12000000)
      throw fpgaHangError("smithwater kernel stuck at redo on fpga");
  }
}

void SWTask::redo_func() {
  ((XCLAgent *)agent_)->fence();
  agent_->writeInput(i_buf[0], i_data[0], i_size[0]*sizeof(int), 0);
  agent_->writeInput(i_buf[1], i_data[1], i_size[1]*sizeof(int), 1);
  agent_->start(this, NULL);
  agent_->readOutput(o_buf[0], o_data[0], o_size[0]*sizeof(int), 0);
  agent_->readOutput(o_buf[1], o_data[1], o_size[1]*sizeof(int), 1);
  agent_->finish();
}

void SWTask::helper_func() {
  while (true) {
    mtx_.lock();
    int exec_state = state_.load();
    if (exec_state == 0) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    }
    else if (exec_state == 1) {
      SWTask *dep_task = this->prv_task_.load();
      this->start_func(dep_task);
      state_.store(0);
    }
    else if (exec_state == 2) {
      this->finish_func();
      state_.store(0);
    }
    else if (exec_state == 3) {
      this->redo_func();
      state_.store(0);
    }
    mtx_.unlock();
  }
}

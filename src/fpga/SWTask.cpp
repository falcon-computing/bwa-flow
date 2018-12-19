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

  i_size = 0;
  o_size = 0;

  i_data = (char*) sw_malloc(max_i_size_, sizeof(char));
  o_data = (short*)sw_malloc(max_o_size_, sizeof(short));

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

#ifndef DEPLOY_aws
  state_.store(0);
  helper_ = boost::thread(boost::bind(&SWTask::helper_func, this));
#endif
}

SWTask::~SWTask() {
#ifndef DEPLOY_aws
  helper_.interrupt();
#endif
  //helper_.detach();
  //pthread_cancel(helper_.native_handle());

  delete region_batch;
  delete chain_batch;

  clReleaseMemObject(i_buf);
  clReleaseMemObject(o_buf);
  free(i_data);
  free(o_data);

  delete agent_;
}

void SWTask::start(SWTask* prev_task) {
  if (o_size == 0) {
    return;
  }

#ifndef DEPLOY_aws
  state_.store(1);
  prv_task_.store(prev_task);

  uint64_t start_ts = getUs();
  while (state_.load() != 0) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (getUs() >= start_ts+10000000) {
      DLOG(ERROR) << "timeout in SWTask::start()";
      throw fpgaHangError("smithwater kernel stuck at start on fpga");
    }
  }
#else
  start_func(prev_task);
#endif
}

void SWTask::start_func(SWTask* prev_task) {
  if (i_size*sizeof(int) >= max_i_size_ || 
      o_size >= max_o_size_) 
  {
    DLOG(ERROR) << "exceeding max memory size";
    throw std::runtime_error("exceeding max memory size");
  }
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Task info: " 
    << "i_size = " << i_size;

  uint64_t start_ts = getUs();
  agent_->writeInput(i_buf, i_data, i_size*sizeof(int), 0);

  if (prev_task->i_size == 0) {
    agent_->start(this, NULL);
  }
  else {
    agent_->start(this, prev_task->agent_);
  }

  agent_->readOutput(o_buf, o_data, o_size*sizeof(int), 0);

  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Enqueue write, kernel & read takes " <<
    getUs() - start_ts << " us";
}

void SWTask::finish() {
  if (o_size == 0) {
    return;
  }

#ifndef DEPLOY_aws
  state_.store(2);

  uint64_t start_ts = getUs();
  while (state_.load() != 0) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (getUs() >= start_ts+10000000) {
      DLOG(ERROR) << "timeout in SWTask::finish()";
      throw fpgaHangError("smithwater kernel stuck at finish on fpga");
    }
  }
#else
  finish_func();
#endif
}

void SWTask::finish_func() {
  uint64_t start_ts = getUs();
  agent_->finish();
  DLOG_IF(INFO, VLOG_IS_ON(4)) << "Dequeue write, kernel & read takes "
                               << getUs() - start_ts << " us";
}

void SWTask::redo() {
  if (o_size == 0) {
    return;
  }

#ifndef DEPLOY_aws
  state_.store(3);

  uint64_t start_ts = getUs();
  while (state_.load() != 0) {
    boost::this_thread::sleep_for(boost::chrono::microseconds(5));
    if (getUs() >= start_ts+10000000) {
      DLOG(ERROR) << "timeout in SWTask::redo()";
      throw fpgaHangError("smithwater kernel stuck at redo on fpga");
    }
  }
#else
  redo_func();
#endif
}

void SWTask::redo_func() {
#ifdef XILINX_FPGA
  ((XCLAgent *)agent_)->fence();
  agent_->writeInput(i_buf, i_data, i_size*sizeof(int), 0);
  agent_->start(this, NULL);
  agent_->readOutput(o_buf, o_data, o_size*sizeof(int), 0);
  agent_->finish();
#endif
}

void sig_handler(int sig) {
  DLOG(INFO) << "caught sig = " << sig;
  throw std::runtime_error("signal caught");
}

void SWTask::helper_func() {
  signal(2, sig_handler);
  while (true) {
    try {
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
    } catch (std::runtime_error & e) {
      DLOG(INFO) << "caught exception";
      break;
    }
  }
}

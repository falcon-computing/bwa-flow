#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>
#include "FPGAAgent.h"

FPGAAgent::FPGAAgent(
    OpenCLEnv* env,
    int chunk_size,
    uint64_t buf_size): 
  env_(env),
  max_buf_size_(buf_size),
  chunk_size_(chunk_size) 
{
  cl_context context = env_->getContext();

  for (int i=0; i<2; i++) {
    // Allocate maximum OpenCL input buffer
    input_buf_[i] = clCreateBuffer(context,
        CL_MEM_READ_ONLY,
        buf_size,
        NULL, NULL);

    // Allocate maximum OpenCL output buffer
    output_buf_[i] = clCreateBuffer(context,
        CL_MEM_WRITE_ONLY, 
        2*chunk_size*FPGA_RET_PARAM_NUM*sizeof(int),
        NULL, NULL);

    fpga_task_[i] = NULL;
  }

  // Skip error handling at this point
}

FPGAAgent::~FPGAAgent() {
  for (int i = 0; i < 2; i++) {
    clReleaseMemObject(input_buf_[i]);
    clReleaseMemObject(output_buf_[i]);
  }
}

void FPGAAgent::writeInput(void* host_ptr, uint64_t size, int cnt) {

  if (size > max_buf_size_) {
    throw std::runtime_error("Input buffer size too big");
  }

  cl_command_queue cmd = env_->getCmdQueue();

  //boost::lock_guard<OpenCLEnv> guard(*env_);

  cl_int err = clEnqueueWriteBuffer(cmd, 
      input_buf_[cnt], CL_TRUE, 0, size, 
      host_ptr, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to write to input buffer");
  }
}

void FPGAAgent::readOutput(void* host_ptr, uint64_t size, int cnt) {
  if (size > 2*chunk_size_*FPGA_RET_PARAM_NUM*sizeof(int)) {
    throw std::runtime_error("Output buffer size too big");
  }
  cl_command_queue cmd = env_->getCmdQueue();
  
  //boost::lock_guard<OpenCLEnv> guard(*env_);

  cl_int err = clEnqueueReadBuffer(cmd,
      output_buf_[cnt], CL_TRUE, 0, size,
      host_ptr, 0, NULL, NULL);  

  if (err != CL_SUCCESS) {
    throw std::runtime_error("Cannot read output buffer\n");
  }
}

void FPGAAgent::start(int task_num, int cnt) {
 // if (task_num > 2*chunk_size_) {
 //   throw std::runtime_error("Too many tasks\n");
 // }

  fpga_task_[cnt] = new FPGATask;
  fpga_task_[cnt]->task_num = task_num;
  fpga_task_[cnt]->input    = &input_buf_[cnt];
  fpga_task_[cnt]->output   = &output_buf_[cnt];

  env_->post_task(fpga_task_[cnt]);
}

void FPGAAgent::wait(int cnt) {
  // Wait for task to finish
  boost::shared_future<bool> f = fpga_task_[cnt]->ready.get_future();

  if (f.valid()) f.wait();

  delete fpga_task_[cnt];
  fpga_task_[cnt] = NULL;
}

bool FPGAAgent::pending(int cnt) {
  if (fpga_task_[cnt]) {
    return true;
  }
  else {
    return false;
  }
}

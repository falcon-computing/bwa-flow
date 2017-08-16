#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "IntelAgent.h"
#include "SWTask.h"

IntelAgent::IntelAgent(OpenCLEnv* env, SWTask* task): task_(task) {

  cl_context context_ = env->getContext();

  cl_int        err       = 0;
  cl_device_id  device_id = env->getDeviceId();
  cl_program    program   = env->getProgram();

  for (int i = 0; i < 4; i++) {
    cmd_[i] = clCreateCommandQueue(context_, device_id, 0, &err);
    if (err != CL_SUCCESS) {
      throw std::runtime_error("Failed to create a command queue context!");
    }
  }
  for (int i = 0; i < 2; i++) {
    char kernel_in_name[100];
    char kernel_out_name[100];

    sprintf(kernel_in_name,"data_parse%d", i);
    sprintf(kernel_out_name,"upload_results%d", i);

    // allocate buffers
    size_t max_i_size = task->max_i_size_;
    size_t max_o_size = task->max_o_size_;;
    i_buf_[i] = clCreateBuffer(context_, CL_MEM_READ_ONLY, sizeof(int)*max_i_size, NULL, NULL);
    o_buf_[i] = clCreateBuffer(context_, CL_MEM_WRITE_ONLY, sizeof(int)*max_o_size, NULL, NULL);

    kernels_[2*i+0] = clCreateKernel(program, kernel_in_name, &err);
    kernels_[2*i+1] = clCreateKernel(program, kernel_out_name, &err);
  }
}

IntelAgent::~IntelAgent() {
  for (int i = 0; i < 2; i++) {
    OCL_CHECKRUN(clReleaseMemObject(i_buf_[i]));
    OCL_CHECKRUN(clReleaseMemObject(o_buf_[i]));
  }
  for (int i = 0; i < 4; i++) {
    OCL_CHECKRUN(clReleaseCommandQueue(cmd_[i]));
    OCL_CHECKRUN(clReleaseKernel(kernels_[i]));
  }
}

void IntelAgent::writeInput(void* host_ptr, int size, int bank) {
  if (size == 0) return;
  uint64_t start_ts = getUs();
  OCL_CHECKRUN(clEnqueueWriteBuffer(cmd_[2*bank], i_buf_[bank],
        CL_FALSE, 0, size, 
        host_ptr, 0, NULL, NULL));
}

void IntelAgent::readOutput(void* host_ptr, int size, int bank) {
  if (size == 0) return;
  OCL_CHECKRUN(clEnqueueReadBuffer(cmd_[2*bank+1], o_buf_[bank],
        CL_TRUE, 0, size, 
        host_ptr, 0, NULL, NULL));
}

void IntelAgent::start(FPGAAgent* agent) {

  IntelAgent* prev_agent = NULL;
  if (agent) {
    prev_agent = (IntelAgent*)agent;
    if (!prev_agent) {
      throw std::runtime_error("unexpected error with FPGA agent");
    }
  }

  // kernel execution
  for (int k = 0; k < 2; k++) {
    cl_int err = 0;
    err  = clSetKernelArg(kernels_[2*k+0], 0, sizeof(cl_mem), &i_buf_[k]);
    err |= clSetKernelArg(kernels_[2*k+0], 1, sizeof(int), &task_->i_size[k]);
    err |= clSetKernelArg(kernels_[2*k+1], 0, sizeof(cl_mem), &o_buf_[k]);
    err |= clSetKernelArg(kernels_[2*k+1], 1, sizeof(int), &task_->o_size[k]);

    if (err) {
      LOG(ERROR) << "failed to set kernel args";
    }
  }
  //cl_event task_event[4];
  for (int k = 0; k < 2; k++) {
    cl_int err = 0;
    if (prev_agent) {
      err  = clEnqueueTask(cmd_[2*k+0], kernels_[2*k+0], 4, prev_agent->events_, &events_[2*k+0]); // input
      err |= clEnqueueTask(cmd_[2*k+1], kernels_[2*k+1], 4, prev_agent->events_, &events_[2*k+1]); // output
    }
    else { // first time the task is executed
      err  = clEnqueueTask(cmd_[2*k+0], kernels_[2*k+0], 0, NULL, &events_[2*k+0]); // input
      err |= clEnqueueTask(cmd_[2*k+1], kernels_[2*k+1], 0, NULL, &events_[2*k+1]); // output
    }

    if (err) {
      LOG(ERROR) << "failed to execute kernels: " << err;
    }
  }
}

void IntelAgent::finish() {
  for (int k = 0; k < 4; k++) {
    clReleaseEvent(events_[k]);
  }
}

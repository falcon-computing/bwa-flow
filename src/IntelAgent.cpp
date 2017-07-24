#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "IntelAgent.h"
#include "SWTask.h"

IntelAgent::IntelAgent(OpenCLEnv* env) {

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

    kernels_[2*i+0] = clCreateKernel(program, kernel_in_name, &err);
    kernels_[2*i+1] = clCreateKernel(program, kernel_out_name, &err);
  }
}

IntelAgent::~IntelAgent() {
  for (int i=0; i<4; i++) {
    clReleaseCommandQueue(cmd_[i]);
    clReleaseKernel(kernels_[i]);
  }
}

void IntelAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  //boost::lock_guard<OpenCLEnv> guard(*env_);
  if (size == 0) return;
  uint64_t start_ts = getUs();
  //cl_event event;
  cl_int err = clEnqueueWriteBuffer(cmd_[2*bank], buf, CL_FALSE, 0, size, 
      host_ptr, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to write buffer to FPGA");
  }
  //err = clFlush(cmd_[2*bank]);
  //if (err != CL_SUCCESS) {
  //  throw std::runtime_error("Failed to flush queue");
  //}
  //clWaitForEvents(1, &event);
  //clReleaseEvent(event);
  //DLOG_IF(INFO, VLOG_IS_ON(3)) << "writeInput:" << bank << 
  //                                ", size=" << size << " takes " <<
  //                                getUs() - start_ts;
}

void IntelAgent::readOutput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size == 0) return;
  cl_int err = clEnqueueReadBuffer(cmd_[2*bank+1], buf, CL_TRUE, 0, size, 
      host_ptr, 0, NULL, NULL);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to read buffer from FPGA");
  }
}

void IntelAgent::start(SWTask* task, FPGAAgent* agent) {

  IntelAgent* prev_agent = NULL;
  if (agent) {
    prev_agent = (IntelAgent*)agent;
    if (!prev_agent) {
      throw std::runtime_error("unexpected error with FPGA agent");
    }
  }

  // kernel execution
  for (int k = 0; k < 2; k++) {
    int i_size = task->i_size[k] / sizeof(int);
    cl_int err = 0;
    err  = clSetKernelArg(kernels_[2*k+0], 0, sizeof(cl_mem), &task->i_buf[k]);
    err |= clSetKernelArg(kernels_[2*k+0], 1, sizeof(int), &i_size);
    err |= clSetKernelArg(kernels_[2*k+1], 0, sizeof(cl_mem), &task->o_buf[k]);
    err |= clSetKernelArg(kernels_[2*k+1], 1, sizeof(int), &task->o_size[k]);

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
  //clWaitForEvents(4, task_event);
  //for (int k = 0; k < 4; k++) {
  //  clReleaseEvent(task_event[k]);
  //}
}

void IntelAgent::finish() {
  for (int k = 0; k < 4; k++) {
    clReleaseEvent(events_[k]);
  }
}

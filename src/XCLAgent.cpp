#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "XCLAgent.h"
#include "SWTask.h"

XCLAgent::XCLAgent(OpenCLEnv* env) {

  cl_context context_ = env->getContext();

  cl_int        err       = 0;
  cl_device_id  device_id = env->getDeviceId();
  cl_program    program   = env->getProgram();

  cmd_ = env->getCmdQueue();
  //cmd_ = clCreateCommandQueue(context_, device_id, 0, &err);
  //if (err != CL_SUCCESS) {
  //  throw std::runtime_error("Failed to create a command queue context");
  //}
  
  kernel_ = clCreateKernel(program, "sw_top", &err);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to create kernel");
  }
}

XCLAgent::~XCLAgent() {
  clReleaseCommandQueue(cmd_);
  clReleaseKernel(kernel_);
}

void XCLAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
    cl_int err = clEnqueueWriteBuffer(cmd_, buf, CL_FALSE, 0, size, 
        host_ptr, 0, NULL, &write_events_[bank]);

    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error writing buffer of size " << size
        << " err: " << err;
      throw std::runtime_error("Failed to write buffer to FPGA");
    }
  }
}

void XCLAgent::readOutput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
    cl_int err = clEnqueueReadBuffer(cmd_, buf, CL_TRUE, 0, size, 
        host_ptr, 1, &kernel_event_, NULL);
    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error reading buffer of size " << size
        << " err: " << err;
      throw std::runtime_error("Failed to read buffer from FPGA");
    }
  }
}

void XCLAgent::start(SWTask* task, FPGAAgent* agent) {

  XCLAgent* prev_agent = NULL;
  if (agent) {
    prev_agent = (XCLAgent*)agent;
    if (!prev_agent) {
      throw std::runtime_error("unexpected error with FPGA agent");
    }
  }

  // kernel execution
  cl_int err = 0;
  int i_arg = 0;
  err  = clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->i_buf[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->i_buf[1]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->o_buf[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->o_buf[1]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &task->i_size[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &task->i_size[1]);

  if (err) {
    LOG(ERROR) << "failed to set kernel args";
  }

  if (prev_agent) {
    cl_event wait_list[3];
    wait_list[0] = prev_agent->kernel_event_;
    wait_list[1] = write_events_[0];
    wait_list[2] = write_events_[1];
    if (task->i_size[1] > 0) {
      err = clEnqueueTask(cmd_, kernel_, 3, wait_list, &kernel_event_);
      valid_2nd_event_ = true;
    }
    else {
      err = clEnqueueTask(cmd_, kernel_, 2, wait_list, &kernel_event_);
      valid_2nd_event_ = false;
    }
  }
  else {
    if (task->i_size[1] > 0) {
      err = clEnqueueTask(cmd_, kernel_, 2, write_events_, &kernel_event_);
      valid_2nd_event_ = true;
    }
    else {
      err = clEnqueueTask(cmd_, kernel_, 1, write_events_, &kernel_event_);
      valid_2nd_event_ = false;
    }
  }
  if (err) {
    LOG(ERROR) << "failed to execute kernels: " << err;
  }
}

void XCLAgent::finish() {
  clReleaseEvent(kernel_event_);
  clReleaseEvent(write_events_[0]);
  if (valid_2nd_event_) {
    clReleaseEvent(write_events_[1]);
  }
}

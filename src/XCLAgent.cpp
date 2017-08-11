#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "XCLAgent.h"
#include "SWTask.h"

XCLAgent::XCLAgent(OpenCLEnv* env, SWTask* task): task_(task) {
  context_ = env->getContext();

  cl_int        err       = 0;
  cl_device_id  device_id = env->getDeviceId();
  cl_program    program   = env->getProgram();

  cmd_ = clCreateCommandQueue(context_, device_id, 0, &err);
  OCL_CHECK(err, "clCreateCommandQueue");
  
  kernel_ = clCreateKernel(program, "sw_top", &err);
  OCL_CHECK(err, "clCreateKernel");
}

XCLAgent::~XCLAgent() {
  clReleaseCommandQueue(cmd_);
  clReleaseKernel(kernel_);
}

void XCLAgent::writeInput(void* host_ptr, int size, int bank) {
  if (size == 0) return;

  cl_int err = 0;
  i_buf_[bank] = clCreateBuffer(context_, 
      CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  
      size, host_ptr, &err);

  OCL_CHECK(err, "clCreateBuffer");

  OCL_CHECKRUN(clEnqueueMigrateMemObjects(cmd_, 1, &i_buf_[bank], 
        0, 0, NULL, &write_events_[bank]));
}

void XCLAgent::readOutput(void* host_ptr, int size, int bank) {
  if (size == 0) return;
  cl_event read_event;
  OCL_CHECKRUN(clEnqueueMigrateMemObjects(cmd_, 1, &o_buf_[bank], 
        CL_MIGRATE_MEM_OBJECT_HOST, 1, &kernel_event_, &read_event));

  OCL_CHECKRUN(clWaitForEvents(1, &read_event));
  OCL_CHECKRUN(clReleaseEvent(read_event));
}

void XCLAgent::start(FPGAAgent* agent) {

  XCLAgent* prev_agent = NULL;
  if (agent) {
    prev_agent = (XCLAgent*)agent;
    if (!prev_agent) {
      throw std::runtime_error("unexpected error with FPGA agent");
    }
  }

  // create output buffer
  for (int k = 0; k < 2; k++) { 
    if (task_->o_size[k] > 0) {
      cl_int err = 0;
      o_buf_[k]  = clCreateBuffer(context_, 
          CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,  
          task_->o_size[k]*sizeof(int), task_->o_data[k], &err);
      OCL_CHECK(err, "clCreateBuffer");
    }
  }

  // kernel execution
  cl_int err = 0;
  int i_arg = 0;
  err  = clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &i_buf_[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &i_buf_[1]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &o_buf_[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &o_buf_[1]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &task_->i_size[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &task_->i_size[1]);

  if (err) {
    LOG(ERROR) << "failed to set kernel args";
  }

  if (prev_agent) {
    cl_event wait_list[3];
    wait_list[0] = prev_agent->kernel_event_;
    wait_list[1] = write_events_[0];
    wait_list[2] = write_events_[1];
    if (task_->i_size[1] > 0) {
      err = clEnqueueTask(cmd_, kernel_, 3, wait_list, &kernel_event_);
    }
    else {
      err = clEnqueueTask(cmd_, kernel_, 2, wait_list, &kernel_event_);
    }
  }
  else {
    if (task_->i_size[1] > 0) {
      err = clEnqueueTask(cmd_, kernel_, 2, write_events_, &kernel_event_);
    }
    else {
      err = clEnqueueTask(cmd_, kernel_, 1, write_events_, &kernel_event_);
    }
  }
  if (err) {
    LOG(ERROR) << "failed to execute kernels: " << err;
  }
}

void XCLAgent::finish() {
  OCL_CHECKRUN(clReleaseEvent(kernel_event_));
  for (int k = 0; k < 2; k++) {
    if (task_->i_size[k]) {
      OCL_CHECKRUN(clReleaseEvent(write_events_[k]));
      OCL_CHECKRUN(clReleaseMemObject(i_buf_[k]));
    }
    if (task_->o_size[k]) OCL_CHECKRUN(clReleaseMemObject(o_buf_[k]));
  }
}

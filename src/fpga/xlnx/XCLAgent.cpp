#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "XCLAgent.h"
#include "SWTask.h"

XCLAgent::XCLAgent(BWAOCLEnv* env, SWTask* task): env_(env) {

  pe_ = env->getPE("sw");

  cl_int     err     = 0;
  cl_context context = pe_.accx->context;
  cl_program program = pe_.accx->program;
  
  kernel_ = clCreateKernel(program, "sw_top", &err);
  OCL_CHECK(err, "failed to create kernel: sw_top");

  cl_mem_ext_ptr_t ext_in, ext_out;

  ext_in.flags = XCL_MEM_DDR_BANK1;
  ext_in.obj = 0; ext_in.param = 0;
  task->i_buf = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX ,
      sizeof(int)*task->max_i_size_, &ext_in, NULL);

  ext_out.flags = XCL_MEM_DDR_BANK1;
  ext_out.obj = 0; ext_out.param = 0;
  task->o_buf = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX ,
      sizeof(int)*task->max_o_size_, &ext_out, NULL);

  kernel_time_ = 0;
  kernel_invks_ = 0;
  writing_time_ = 0;
  reading_time_ = 0;
}

XCLAgent::~XCLAgent() {
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Kernel total time: " << kernel_time_ << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Kernel average time: " << (double)kernel_time_/(double)kernel_invks_ << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Writing total time: " << writing_time_ << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[sw] Reading total time: " << reading_time_ << " us";
  clReleaseKernel(kernel_);
  env_->releasePE(pe_);
}

void XCLAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  //cl_command_queue cmd = pe_.accx->cmd;
  cl_command_queue cmd = pe_.cmd;
  if (size > 0) {
    cl_int err = clEnqueueWriteBuffer(cmd, buf, CL_FALSE, 0, size, host_ptr, 0, NULL, &write_events_);
    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error writing buffer of size " << size
        << " err: " << err;
      throw std::runtime_error("Failed to write buffer to FPGA");
    }
  }
}

void XCLAgent::readOutput(cl_mem buf, void* host_ptr, int size, int bank) {
  //cl_command_queue cmd = pe_.accx->cmd;
  cl_command_queue cmd = pe_.cmd;
  if (size > 0) {
    cl_int err = clEnqueueReadBuffer(cmd, buf, CL_FALSE, 0, size, host_ptr, 1, &kernel_event_, &read_events_);
    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error reading buffer of size " << size
        << " err: " << err;
      throw std::runtime_error("Failed to read buffer from FPGA");
    }
  }
}

void XCLAgent::start(Task* i_task, FPGAAgent* agent) {

  SWTask *task = (SWTask *)i_task;
  cl_command_queue cmd = pe_.cmd;

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
  err  = clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->i_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &task->o_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &env_->pac_input_list_[pe_.accx->accx_id]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &task->i_size);

  if (err) {
    LOG(ERROR) << "failed to set kernel args";
  }

  if (prev_agent) {
    cl_event wait_list[2];
    wait_list[0] = prev_agent->kernel_event_;
    wait_list[1] = write_events_;
    err = clEnqueueTask(cmd, kernel_, 2, wait_list, &kernel_event_);
  }
  else {
    err = clEnqueueTask(cmd, kernel_, 1, &write_events_, &kernel_event_);
  }
  if (err) {
    LOG(ERROR) << "failed to execute kernels: " << err;
  }
  //clFlush(cmd);
}

void XCLAgent::finish() {
  clWaitForEvents(1, &read_events_);

#if 0
  cl_ulong k_start, k_end;
  clGetEventProfilingInfo(write_events_[0], CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
  clGetEventProfilingInfo(write_events_[0], CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
  writing_time_ += (uint64_t)(k_end - k_start)/1000;
  if (valid_2nd_event_) {
    clGetEventProfilingInfo(write_events_[1], CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
    clGetEventProfilingInfo(write_events_[1], CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
    writing_time_ += (uint64_t)(k_end - k_start)/1000;
  }

  clGetEventProfilingInfo(kernel_event_, CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
  clGetEventProfilingInfo(kernel_event_, CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
  kernel_time_ += (uint64_t)(k_end - k_start)/1000;
  kernel_invks_++;

  clGetEventProfilingInfo(read_events_[0], CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
  clGetEventProfilingInfo(read_events_[0], CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
  reading_time_ += (uint64_t)(k_end - k_start)/1000;
  if (valid_2nd_event_) {
    clGetEventProfilingInfo(read_events_[1], CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
    clGetEventProfilingInfo(read_events_[1], CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
    reading_time_ += (uint64_t)(k_end - k_start)/1000;
  }
#endif

  clReleaseEvent(write_events_);
  clReleaseEvent(kernel_event_);
  clReleaseEvent(read_events_);
}

void XCLAgent::fence() {
  //cl_command_queue cmd = pe_.accx->cmd;
  cl_command_queue cmd = pe_.cmd;
  clFlush(cmd);
  clFinish(cmd);
}

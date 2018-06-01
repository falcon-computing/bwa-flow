#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "SMemXCLAgent.h"
#include "SMemTask.h"


SMemXCLAgent::SMemXCLAgent(BWAOCLEnv* env, SMemTask* task) {
  env_ = env;

  cl_int err = 0;
  cl_program program = env->getProgram();
  kernel_ = clCreateKernel(program, "mem_collect_intv_fpga", &err);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to create kernel");
  }
}

SMemXCLAgent::~SMemXCLAgent() {
  clReleaseKernel(kernel_);
}

void SMemXCLAgent::createBuffer(SMemTask *task) {
  cl_context context = env_->getContext();
  task->i_seq_buf = clCreateBuffer( context, 
                                    CL_MEM_READ_ONLY, 
                                    sizeof(uint8_t)*task->max_i_seq_len_*task->max_i_seq_num_, 
                                    NULL, 
                                    NULL );
  cl_mem_ext_ptr_t ext;
  ext.flags = XCL_MEM_DDR_BANK0; ext.obj = 0; ext.param=0;
  task->o_mem_buf = clCreateBuffer( context,
                                    CL_MEM_WRITE_ONLY | CL_MEM_EXT_PTR_XILINX, 
                                    sizeof(bwtintv_t)*task->max_intv_alloc_*task->max_i_seq_num_,
                                    &ext,
                                    NULL );
  task->o_num_buf = clCreateBuffer( context, 
                                    CL_MEM_WRITE_ONLY|CL_MEM_EXT_PTR_XILINX, 
                                    sizeof(int)*task->max_i_seq_num_,
                                    &ext,
                                    NULL );
} 

void SMemXCLAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
    cl_command_queue command = env_->getCmdQueue();
    cl_int err = clEnqueueWriteBuffer( command, buf, CL_FALSE, 0,
                                       size, host_ptr, 
                                       0, NULL, &write_events_ );
    
    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error writing buffer of size " << size
                  << " err: " << err;
      throw std::runtime_error("Failed to write buffer to FPGA");
    }
  }
}

void SMemXCLAgent::readOutput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
#if 0
    cl_int status;
    cl_int e;
    e = clGetEventInfo( kernel_event_, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &status, NULL );
    if      ( e != CL_SUCCESS        )  DLOG(INFO) << "Invalid.";
    else if ( status == CL_QUEUED    )  DLOG(INFO) << "Queued.";
    else if ( status == CL_SUBMITTED )  DLOG(INFO) << "Submitted.";
    else if ( status == CL_RUNNING   )  DLOG(INFO) << "Running.";
    else if ( status == CL_COMPLETE  )  DLOG(INFO) << "Complete.";
    else                                DLOG(INFO) << "Unknown status";
#endif
    cl_command_queue command = env_->getCmdQueue();
    cl_int err = clEnqueueReadBuffer( command, buf, CL_TRUE, 0,
                                      size, host_ptr,
                                      1, &kernel_event_, NULL );
    if (err != CL_SUCCESS) {
      DLOG(ERROR) << "error reading buffer of size " << size
                  << " err: " << err;
      throw std::runtime_error("Failed to read buffer from FPGA");
    }
  }
}

void SMemXCLAgent::start(Task* task, FPGAAgent* agent) {
  SMemTask *cur_task = (SMemTask *)task;
  SMemXCLAgent* prev_agent = NULL;
  if (agent) {
    prev_agent = (SMemXCLAgent *)agent;
    if (!prev_agent) {
      throw std::runtime_error("unexpected error with FPGA agent");
    }
  }

  // set kernel args
  cl_int err = 0;
  int i_arg = 0;
  err  = clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &env_->bwt_[0]);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->i_seq_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->o_mem_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->o_num_buf);
  for (int i=1; i<16; i++) {
    err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &env_->bwt_[i]);
  }
  if (err) {
    LOG(ERROR) << "failed to set kernel args";
    throw std::runtime_error("Failed to set kernel args");
  }

  // enqueue kernel
  cl_command_queue command = env_->getCmdQueue();
  if (prev_agent) {
    cl_event wait_list[2];
    wait_list[0] = prev_agent->kernel_event_;
    wait_list[1] = write_events_;
    err = clEnqueueTask(command, kernel_, 2, wait_list, &kernel_event_);
    
  }
  else {
    err = clEnqueueTask(command, kernel_, 1, &write_events_, &kernel_event_);
  }
  if (err) {
    LOG(ERROR) << "failed to execute kernels: " << err;
    throw std::runtime_error("Failed to execute the kernel");
  }
}

void SMemXCLAgent::finish() {
  clReleaseEvent(kernel_event_);
  clReleaseEvent(write_events_);
}

void SMemXCLAgent::releaseBuffer(SMemTask *task) {
  clReleaseMemObject(task->i_seq_buf);
  clReleaseMemObject(task->o_mem_buf);
  clReleaseMemObject(task->o_num_buf);
}

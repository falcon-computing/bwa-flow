#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "SMemXCLAgent.h"
#include "SMemTask.h"


SMemXCLAgent::SMemXCLAgent(BWAOCLEnv* env, SMemTask* task) {
  env_ = env;
  device_ = env->getDevice();

  cl_int err = 0;
  cl_context context = device_.context;
  cl_program program = device_.program;

  for (int bank_id=0; bank_id<SMEM_BANK_NUM; bank_id++)
    kernel_[bank_id] = clCreateKernel(program, std::string("mem_collect_intv_core"+std::to_string(bank_id)).c_str(), &err);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to create kernel");
  }
}

SMemXCLAgent::~SMemXCLAgent() {
  for (int bank_id=0; bank_id<SMEM_BANK_NUM; bank_id++)
    clReleaseKernel(kernel_[bank_id]);
  env_->releaseDevice(device_);
}

void SMemXCLAgent::createBuffer(SMemTask *task) {
  static const unsigned bankID[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
  cl_context context = device_.context;
  for (int i=0; i<SMEM_BANK_NUM; i++) {
    cl_mem_ext_ptr_t ext;
    ext.flags = bankID[i];
    ext.obj = 0; ext.param=0;
    task->i_seq_buf[i] = clCreateBuffer( context, 
                                         CL_MEM_READ_ONLY|CL_MEM_EXT_PTR_XILINX, 
                                         task->i_seq_size[i],/*sizeof(uint8_t)*task->max_i_seq_len_*task->max_i_seq_num_, */
                                         &ext, 
                                         NULL );
    task->i_seq_len_buf[i] = clCreateBuffer( context, 
                                             CL_MEM_READ_ONLY|CL_MEM_EXT_PTR_XILINX, 
                                             task->i_seq_len_size[i],/*sizeof(uint8_t)*task->max_i_seq_num_, */
                                             &ext, 
                                             NULL );
    task->o_mem_buf[i] = clCreateBuffer( context,
                                         CL_MEM_WRITE_ONLY|CL_MEM_EXT_PTR_XILINX, 
                                         task->o_mem_size[i],/*sizeof(bwtintv_t)*task->max_intv_alloc_*task->max_i_seq_num_,*/
                                         &ext,
                                         NULL );
    task->o_num_buf[i] = clCreateBuffer( context, 
                                         CL_MEM_WRITE_ONLY|CL_MEM_EXT_PTR_XILINX, 
                                         task->o_num_size[i],/*sizeof(int)*task->max_i_seq_num_,*/
                                         &ext,
                                         NULL );
  }
} 

void SMemXCLAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
    cl_command_queue cmd = device_.cmd;
    // non-blocking writing input
    //cl_int err = clEnqueueWriteBuffer( cmd, buf, CL_FALSE, 0,
    //                                   size, host_ptr, 
    //                                   0, NULL, &write_events_[bank] );

    // blocking writing input
    cl_int err = clEnqueueWriteBuffer( cmd, buf, CL_TRUE, 0,
                                       size, host_ptr, 
                                       0, NULL, NULL );

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
    cl_command_queue cmd = device_.cmd;
    cl_int err = clEnqueueReadBuffer( cmd, buf, CL_TRUE, 0,
                                      size, host_ptr,
                                      1, &kernel_event_[bank], NULL );
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

  for (int i = 0; i < SMEM_BANK_NUM; i++) {
    // set kernel args
    cl_int err = 0;
    int i_arg = 0;
    err  = clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &env_->bwt_list_[device_.env_id][i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &env_->bwt_param_list_[device_.env_id][i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &cur_task->i_seq_buf[i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &cur_task->i_seq_len_buf[i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &cur_task->o_mem_buf[i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(cl_mem), &cur_task->o_num_buf[i]);
    err |= clSetKernelArg(kernel_[i], i_arg++, sizeof(int), &cur_task->i_seq_num[i]);
    if (err) {
      LOG(ERROR) << "failed to set kernel args";
      throw std::runtime_error("Failed to set kernel args");
    }

    // enqueue kernel
    cl_command_queue cmd = device_.cmd;
    if (prev_agent) {
      // non-blocking writing input
      //cl_event wait_list[2];
      //wait_list[0] = prev_agent->kernel_event_[i];
      //wait_list[1] = write_events_[i];
      //err = clEnqueueTask(cmd, kernel_[i], 2, wait_list, &kernel_event_[i]);
      
      // blocking writing input
      err = clEnqueueTask(cmd, kernel_[i], 1, &(prev_agent->kernel_event_[i]), &kernel_event_[i]);
      
    }
    else {
      // non-blocking writing input
      //err = clEnqueueTask(cmd, kernel_[i], 1, &write_events_[i], &kernel_event_[i]);

      // blocking writing input
      err = clEnqueueTask(cmd, kernel_[i], 0, NULL, &kernel_event_[i]);
    }
    if (err) {
      LOG(ERROR) << "failed to execute kernels: " << err;
      throw std::runtime_error("Failed to execute the kernel "+std::to_string(i));
    }
  }
}

void SMemXCLAgent::finish() {
  for (int i=0; i<SMEM_BANK_NUM; i++) {
    clReleaseEvent(kernel_event_[i]);
    clReleaseEvent(write_events_[i]);
  }
}

void SMemXCLAgent::releaseBuffer(SMemTask *task) {
  for (int i=0; i<SMEM_BANK_NUM; i++) {
    clReleaseMemObject(task->i_seq_buf[i]);
    clReleaseMemObject(task->i_seq_len_buf[i]);
    clReleaseMemObject(task->o_mem_buf[i]);
    clReleaseMemObject(task->o_num_buf[i]);
  }
}

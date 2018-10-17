#include <boost/thread/future.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

#include "SMemXCLAgent.h"
#include "SMemTask.h"



SMemXCLAgent::SMemXCLAgent(BWAOCLEnv* env, SMemTask* task) {
  const unsigned ddr_bank[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
  env_ = env;
  pe_ = env->getPE("smem");

  cl_int err = 0;
  cl_context context = pe_.accx->context;
  cl_program program = pe_.accx->program;

  kernel_ = clCreateKernel(program, std::string("mem_collect_intv_core"+std::to_string(pe_.bank_id)).c_str(), &err);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to create kernel");
  }

  task->i_bwt = env->bwt_list_[pe_.pe_id];
  task->i_bwt_param = env->bwt_param_list_[pe_.pe_id]; 
  DLOG(INFO) << "pe id: " << pe_.pe_id << "; bank id: " << pe_.bank_id;
  kernel_time_ = 0;
  kernel_invks_ = 0;
}

SMemXCLAgent::~SMemXCLAgent() {
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Kernel total time: " << kernel_time_ << " us";
  DLOG_IF(INFO, VLOG_IS_ON(2)) << "[smem] Kernel average time: " << (double)kernel_time_/(double)kernel_invks_ << " us";
  clReleaseKernel(kernel_);
  env_->releasePE(pe_);
}

void SMemXCLAgent::createBuffer(SMemTask *task) {
  const unsigned ddr_bank[4] = {XCL_MEM_DDR_BANK0, XCL_MEM_DDR_BANK1, XCL_MEM_DDR_BANK2, XCL_MEM_DDR_BANK3};
  cl_context context = pe_.accx->context;
  cl_mem_ext_ptr_t ext_in_seq, ext_in_len, ext_out_mem, ext_out_num;
  ext_in_seq.flags = ddr_bank[pe_.bank_id]; ext_in_seq.obj = task->i_seq_data; ext_in_seq.param=0;
  task->i_seq_buf = clCreateBuffer( context, 
                                    CL_MEM_READ_ONLY|CL_MEM_EXT_PTR_XILINX|CL_MEM_USE_HOST_PTR, 
                                    task->i_seq_size,/*sizeof(uint8_t)*task->max_i_seq_len_*task->max_i_seq_num_, */
                                    &ext_in_seq, 
                                    NULL );
  ext_in_len.flags = ddr_bank[pe_.bank_id]; ext_in_len.obj = task->i_seq_len_data; ext_in_len.param=0;
  task->i_seq_len_buf = clCreateBuffer( context, 
                                        CL_MEM_READ_ONLY|CL_MEM_EXT_PTR_XILINX|CL_MEM_USE_HOST_PTR, 
                                        task->i_seq_len_size,/*sizeof(uint8_t)*task->max_i_seq_num_, */
                                        &ext_in_len, 
                                        NULL );
  ext_out_mem.flags = ddr_bank[pe_.bank_id]; ext_out_mem.obj = task->o_mem_data; ext_out_mem.param=0;
  task->o_mem_buf = clCreateBuffer( context,
                                    CL_MEM_WRITE_ONLY|CL_MEM_EXT_PTR_XILINX|CL_MEM_USE_HOST_PTR, 
                                    task->o_mem_size,/*sizeof(bwtintv_t)*task->max_intv_alloc_*task->max_i_seq_num_,*/
                                    &ext_out_mem,
                                    NULL );
  ext_out_num.flags = ddr_bank[pe_.bank_id]; ext_out_num.obj = task->o_num_data; ext_out_num.param=0;
  task->o_num_buf = clCreateBuffer( context, 
                                    CL_MEM_WRITE_ONLY|CL_MEM_EXT_PTR_XILINX|CL_MEM_USE_HOST_PTR, 
                                    task->o_num_size,/*sizeof(int)*task->max_i_seq_num_,*/
                                    &ext_out_num,
                                    NULL );
} 

void SMemXCLAgent::writeInput(cl_mem buf, void* host_ptr, int size, int bank) {
  if (size > 0) {
    //cl_command_queue cmd = pe_.accx->cmd;
    cl_command_queue cmd = pe_.cmd;
    // non-blocking writing input
    cl_int err = clEnqueueWriteBuffer( cmd, buf, CL_FALSE, 0,
                                       size, host_ptr, 
                                       0, NULL, &write_events_[bank] );
    //cl_int err = clEnqueueMigrateMemObjects(cmd, 1, &buf, 0, 0, NULL, &write_events_[bank] );

    //// blocking writing input
    //cl_int err = clEnqueueWriteBuffer( cmd, buf, CL_TRUE, 0,
    //                                   size, host_ptr, 
    //                                   0, NULL, NULL );

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
    //cl_command_queue cmd = pe_.accx->cmd;
    cl_command_queue cmd = pe_.cmd;
    cl_int err = clEnqueueReadBuffer( cmd, buf, CL_TRUE, 0,
                                      size, host_ptr,
                                      1, &kernel_event_, NULL );
    //cl_event read;
    //cl_int err = clEnqueueMigrateMemObjects(cmd, 1, &buf, CL_MIGRATE_MEM_OBJECT_HOST, 1, &kernel_event_, &read);
    //clWaitForEvents(1, &read);
    //clReleaseEvent(read);
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
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->i_bwt_param);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->i_seq_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->i_seq_len_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->o_mem_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->o_num_buf);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(int), &cur_task->i_seq_num);
  err |= clSetKernelArg(kernel_, i_arg++, sizeof(cl_mem), &cur_task->i_bwt);
  if (err) {
    LOG(ERROR) << "failed to set kernel args";
    throw std::runtime_error("Failed to set kernel args");
  }

  // enqueue kernel
  //cl_command_queue cmd = pe_.accx->cmd;
  cl_command_queue cmd = pe_.cmd;
  if (prev_agent) {
    //DLOG(INFO) << "with prev";
    // non-blocking writing input
    cl_event wait_list[3];
    wait_list[0] = prev_agent->kernel_event_;
    wait_list[1] = write_events_[0];
    wait_list[2] = write_events_[1];
    err = clEnqueueTask(cmd, kernel_, 3, wait_list, &kernel_event_);
    
    //// blocking writing input
    //err = clEnqueueTask(cmd, kernel_, 1, &(prev_agent->kernel_event_), &kernel_event_);
    
  }
  else {
    //DLOG(INFO) << "no prev";
    // non-blocking writing input
    err = clEnqueueTask(cmd, kernel_, 2, write_events_, &kernel_event_);

    //// blocking writing input
    //err = clEnqueueTask(cmd, kernel_, 0, NULL, &kernel_event_);
  }
  if (err) {
    LOG(ERROR) << "failed to execute kernels: " << err;
    throw std::runtime_error("Failed to execute the kernel "+std::to_string(pe_.bank_id));
  }
  
}

void SMemXCLAgent::finish() {
  clReleaseEvent(kernel_event_);
  clReleaseEvent(write_events_[0]);
  clReleaseEvent(write_events_[1]);
}

void SMemXCLAgent::releaseBuffer(SMemTask *task) {
  clReleaseMemObject(task->i_seq_buf);
  clReleaseMemObject(task->i_seq_len_buf);
  clReleaseMemObject(task->o_mem_buf);
  clReleaseMemObject(task->o_num_buf);
}

void SMemXCLAgent::wait() {
  clWaitForEvents(1, &kernel_event_);
#if 0
  cl_ulong k_start, k_end;
  clGetEventProfilingInfo(kernel_event_, CL_PROFILING_COMMAND_START, sizeof(k_start), &k_start, NULL);
  clGetEventProfilingInfo(kernel_event_, CL_PROFILING_COMMAND_END, sizeof(k_end), &k_end, NULL);
  kernel_time_ += (uint64_t)(k_end - k_start)/1000;
  kernel_invks_++;
#endif
}

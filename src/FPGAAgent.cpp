#include "FPGAAgent.h"

FPGAAgent::FPGAAgent(
    const char* bit_path,
    int chunk_size,
    uint64_t buf_size): 
  env_(bit_path, "sw_top"),
  max_buf_size_(buf_size),
  chunk_size_(chunk_size) 
{
  cl_context context = env_.getContext();

  for (int i=0; i<2; i++) {
    // Allocate maximum OpenCL input buffer
    input_buf_[i] = clCreateBuffer(context,
        CL_MEM_READ_ONLY,
        buf_size,
        NULL, NULL);

    // Allocate maximum OpenCL output buffer
    output_buf_[i] = clCreateBuffer(context,
        CL_MEM_WRITE_ONLY, 
        chunk_size*FPGA_RET_PARAM_NUM*sizeof(int),
        NULL, NULL);
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

  cl_command_queue cmd = env_.getCmdQueue();
  cl_int err = clEnqueueWriteBuffer(cmd, 
      input_buf_[cnt], CL_TRUE, 0, size, 
      host_ptr, 0, NULL, NULL);

  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to write to input buffer");
  }
}

void FPGAAgent::readOutput(void* host_ptr, uint64_t size, int cnt) {
  if (size > chunk_size_*FPGA_RET_PARAM_NUM*sizeof(int)) {
    throw std::runtime_error("Output buffer size too big");
  }
  cl_command_queue cmd = env_.getCmdQueue();
  
  cl_int err = clEnqueueReadBuffer(cmd,
      output_buf_[cnt], CL_TRUE, 0, size,
      host_ptr, 0, NULL, NULL);  

  if (err != CL_SUCCESS) {
    throw std::runtime_error("Cannot read output buffer\n");
  }
}

void FPGAAgent::start(int task_num, int cnt) {
  cl_event event;
  cl_int err;
  cl_kernel kernel = env_.getKernel();
  cl_command_queue cmd = env_.getCmdQueue();

  if (task_num > chunk_size_) {
    throw std::runtime_error("Too many tasks\n");
  }

  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_buf_[cnt]);
  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buf_[cnt]);
  err |= clSetKernelArg(kernel, 2, sizeof(int), &task_num);
  if (err != CL_SUCCESS) {
    throw std::runtime_error("Failed to set kernel args!");
  }

  err = clEnqueueTask(cmd, kernel, 0, NULL, &event);
  if (err) {
    throw("Failed to execute kernel!");
  }
  clWaitForEvents(1, &event);
}

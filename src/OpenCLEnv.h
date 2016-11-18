#ifndef OPENCLENV_H
#define OPENCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <glog/logging.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

#include "kflow/Queue.h"
#include "util.h"

struct FPGATask {
  int     size_a;
  int     size_b;
  cl_mem* input_a;
  cl_mem* input_b;
  cl_mem* output_a;
  cl_mem* output_b;
  boost::promise<bool> ready;
};

class OpenCLEnv 
: public boost::basic_lockable_adapter<boost::mutex> {
 public:
  OpenCLEnv(const char* bin_path, const char* kernel_name) {
    // start platform setting up
    int err;

    cl_platform_id* platform_id = new cl_platform_id[2];
    cl_device_id device_id;

    char cl_platform_vendor[1001];
    char cl_platform_name[1001];

    cl_uint num_platforms = 0;

    // Connect to first platform
    err = clGetPlatformIDs(2, platform_id, &num_platforms);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "Failed to find an OpenCL platform!");
    }

    err = clGetPlatformInfo(
        platform_id[1], 
        CL_PLATFORM_VENDOR, 
        1000, 
        (void *)cl_platform_vendor,NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error(
            "clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!");
    }
    err = clGetPlatformInfo(platform_id[1],CL_PLATFORM_NAME,1000,(void *)cl_platform_name,NULL);
    
    if (err != CL_SUCCESS) {
        throw std::runtime_error("clGetPlatformInfo(CL_PLATFORM_NAME) failed!");
    }

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id[1], CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group!");
    }

    // Create a compute context 
    context_ = clCreateContext(0, 1, &device_id, NULL, NULL, &err);

    if (!context_) {
        throw std::runtime_error("Failed to create a compute context!");
    }

    // Create a command commands
    commands_ = clCreateCommandQueue(context_, device_id, 0, &err);
    //commands_ = clCreateCommandQueue(context_, device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);

    if (!commands_) {
        throw std::runtime_error("Failed to create a command queue context!");
    }

    // Create Program Objects
    // TODO: this part should not be static, the program
    // should be configurable at runtime

    // Load binary from disk
    unsigned char *kernelbinary;
    
    int n_i = load_file(bin_path, (char **) &kernelbinary);

    if (n_i < 0) {
      throw std::runtime_error(
          "failed to load kernel from xclbin");
    }
    size_t n_t = n_i;

    int status = 0;

    // Create the compute program from offline
    program_ = clCreateProgramWithBinary(context_, 1, &device_id, &n_t,
            (const unsigned char **) &kernelbinary, &status, &err);

    if ((!program_) || (err!=CL_SUCCESS)) {
        throw std::runtime_error(
            "Failed to create compute program from binary");
    }

    // Build the program executable
    err = clBuildProgram(program_, 0, NULL, NULL, NULL, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to build program executable!");
    }

    // Create the compute kernel in the program we wish to run
    kernel_ = clCreateKernel(program_, kernel_name, &err);

    if (!kernel_ || err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create compute kernel!");
    }

    // Start executor
    task_workers_.create_thread(boost::bind(&OpenCLEnv::execute, this));
  }

  ~OpenCLEnv() {
    task_workers_.interrupt_all();
    task_workers_.join_all();

    clReleaseKernel(kernel_);
    clReleaseCommandQueue(commands_);
    clReleaseProgram(program_);
    clReleaseContext(context_);
  }

  cl_context& getContext() {
    return context_;
  }

  cl_command_queue& getCmdQueue() {
    return commands_;
  }

  cl_kernel& getKernel() {
    return kernel_;
  }

  void post_task(FPGATask* task) {
    task_queue_.push(task);   
  }

  void execute() {
    DLOG(INFO) << "OpenCLEnv executor started.";

    cl_int err;
    cl_event event;

    int counter = 0;
    int thresh_counter = 0;
    uint64_t total_wait_time = 0;
    uint64_t total_fpga_time = 0;

    try {
      while (true) {
        uint64_t start_ts = getNs();

        FPGATask* task;
        task_queue_.pop(task);
        uint64_t wait_time = getNs() - start_ts;

        start_ts = getNs();

        // Got new task and start execution
        int size_a = task->size_a;
        int size_b = task->size_b;
        cl_mem* input_a  = task->input_a;
        cl_mem* output_a = task->output_a;
        cl_mem* input_b  = task->input_b;
        cl_mem* output_b = task->output_b;
        // to make sure the fpga works
        if (size_b == 0) {
          input_b = input_a;
        }
          
        err  = clSetKernelArg(kernel_, 0, sizeof(cl_mem), input_a);
        err |= clSetKernelArg(kernel_, 1, sizeof(cl_mem), input_b);
        err |= clSetKernelArg(kernel_, 2, sizeof(cl_mem), output_a);
        err |= clSetKernelArg(kernel_, 3, sizeof(cl_mem), output_b);
        err |= clSetKernelArg(kernel_, 4, sizeof(int), &size_a);
        err |= clSetKernelArg(kernel_, 5, sizeof(int), &size_b);
        err = clEnqueueTask(commands_, kernel_, 0, NULL, &event);
        if (err) {
          LOG(ERROR) << "Failed to execute kernel.";
        }
        clWaitForEvents(1, &event);

        task->ready.set_value(true);

        uint64_t fpga_time = getNs() - start_ts;

        VLOG(3) << "SW-FPGA kernel takes " 
          << fpga_time/1e3 << " us";

        // Only measure middle part of execution
        if (counter > 100 && counter < 10000) {
          total_wait_time += wait_time;
          total_fpga_time += fpga_time;
        }

        if (wait_time < 0.2*fpga_time) {
          thresh_counter++;
        }
        counter++;
      }
    }
    catch (boost::thread_interrupted &e) {
      VLOG(1) << "FPGA utilization is " 
        << 100.0f*total_fpga_time/(total_fpga_time+total_wait_time) << " %";
      VLOG(1) << "Percentage of tasks that FPGA util > 80\% is " 
        << 100.0f*thresh_counter/counter << " %";
    }
  }

 private:
  int load_file(
      const char *filename, 
      char **result) { 
    int size = 0;
    FILE *f = fopen(filename, "rb");
    if (f == NULL) 
    { 
      *result = NULL;
      return -1; // -1 means file opening fail 
    } 
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) 
    { 
      free(*result);
      return -2; // -2 means file reading fail 
    } 
    fclose(f);
    (*result)[size] = 0;
    return size;
  }

  boost::thread_group task_workers_;
  kestrelFlow::Queue<FPGATask*, 4> task_queue_;

  cl_context       context_;   // compute context
  cl_command_queue commands_;  // compute command queue
  cl_program       program_;   // compute program
  cl_kernel        kernel_;    // compute kernel
};
#endif

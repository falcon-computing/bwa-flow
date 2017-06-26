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
  int     out_size_a;
  int     out_size_b;
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
    DLOG(INFO) << "Found " << num_platforms << " opencl platforms";

    int platform_idx;
    for (int i = 0; i < num_platforms; i++) {
        char cl_platform_name[1001];

        err = clGetPlatformInfo(
            platform_id[i], 
            CL_PLATFORM_NAME, 
            1000, 
            (void *)cl_platform_name,NULL);

        DLOG(INFO) << "Found platform " << cl_platform_name;

        if (err != CL_SUCCESS) {
          throw std::runtime_error(
              "clGetPlatformInfo(CL_PLATFORM_NAME) failed!");
        }

        if (strstr(cl_platform_name, "Xilinx") != NULL || strstr(cl_platform_name, "Intel") != NULL) {
            // found platform
            //printf("Found Xilinx platform\n");
            platform_idx = i;
            DLOG(INFO) << "Selecting platform " << cl_platform_name;
            break;
        }
    }

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR, 1, &device_id, NULL);

    if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group!");
    }

    // Create a compute context 
    context_ = clCreateContext(0, 1, &device_id, NULL, NULL, &err);

    if (!context_) {
      throw std::runtime_error("Failed to create a compute context!");
    }

    // Create a command queues
    for (int i = 0; i < 8; i++) {
        commands_[i] = clCreateCommandQueue(context_, device_id, 0, &err);

        if (err != CL_SUCCESS) {
          throw std::runtime_error("Failed to create a command queue context!");
        }
    }

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
    for (int i = 0; i < 2; i++) {
      char kernel_in_name[100];
      char kernel_out_name[100];

      sprintf(kernel_in_name,"data_parse%d", i);
      sprintf(kernel_out_name,"upload_results%d", i);

      kernels_[2*i+0] = clCreateKernel(program_, kernel_in_name, &err);
      kernels_[2*i+1] = clCreateKernel(program_, kernel_out_name, &err);
    }

    // Start executor
    task_workers_.create_thread(boost::bind(&OpenCLEnv::execute, this));
  }

  ~OpenCLEnv() {
    task_workers_.interrupt_all();
    task_workers_.join_all();

    for (int i = 0; i < 4; i++) {
      clReleaseKernel(kernels_[i]);
      clReleaseCommandQueue(commands_[i]);
    }
    for (int i = 4; i < 8; i++) {
      clReleaseCommandQueue(commands_[i]);
    }
    clReleaseProgram(program_);
    clReleaseContext(context_);
  }

  cl_context& getContext() {
    return context_;
  }

  cl_command_queue& getCmdQueue(int idx) {
    return commands_[idx];
  }

  //cl_kernel& getKernel() {
  //  return kernel_;
  //}

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

    // setup kernels
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
        int out_size_a = task->out_size_a;
        int out_size_b = task->out_size_b;
        cl_mem* input_a  = task->input_a;
        cl_mem* output_a = task->output_a;
        cl_mem* input_b  = task->input_b;
        cl_mem* output_b = task->output_b;
        // to make sure the fpga works
        if (size_b == 0) {
          input_b = input_a;
        }
          
        // kernel execution
        err  = clSetKernelArg(kernels_[0], 0, sizeof(cl_mem), input_a);
        err |= clSetKernelArg(kernels_[0], 1, sizeof(int), &size_a);
        err |= clSetKernelArg(kernels_[2], 0, sizeof(cl_mem), input_b);
        err |= clSetKernelArg(kernels_[2], 1, sizeof(int), &size_b);

        err |= clSetKernelArg(kernels_[1], 0, sizeof(cl_mem), output_a);
        err |= clSetKernelArg(kernels_[1], 1, sizeof(int), &out_size_a);
        err |= clSetKernelArg(kernels_[3], 0, sizeof(cl_mem), output_b);
        err |= clSetKernelArg(kernels_[3], 1, sizeof(int), &out_size_b);

        //boost::lock_guard<OpenCLEnv> guard(this);
        cl_event kernel_events[4];
        for (int i = 0; i < 4; i++) {
          err = clEnqueueTask(commands_[i], kernels_[i], 0, NULL, 
              &kernel_events[i]);

          if (err) {
            LOG(ERROR) << "Failed to execute kernel " << i;
          }
        }
        //clWaitForEvents(1, &event);
        err = clWaitForEvents(4, kernel_events);
        if (err != CL_SUCCESS) {
          LOG(ERROR) << "Failed to wait for events.";
        }

        // release cl_event
        for (int i = 0; i < 4; i++) {
          clReleaseEvent(kernel_events[i]);
        }

        task->ready.set_value(true);

        uint64_t fpga_time = getNs() - start_ts;

        DLOG_IF(INFO, FLAGS_v >= 3) << "SW-FPGA kernel takes " 
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
      DLOG_IF(INFO, FLAGS_v >= 1) << "FPGA utilization is " 
        << 100.0f*total_fpga_time/(total_fpga_time+total_wait_time) << " %";
      DLOG_IF(INFO, FLAGS_v >= 1) << "Percentage of tasks that FPGA util > 80\% is " 
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

  cl_context       context_;      // compute context
  cl_command_queue commands_[8];  // compute command queue
  cl_program       program_;      // compute program
  cl_kernel        kernels_[4];   // compute kernel
};
#endif

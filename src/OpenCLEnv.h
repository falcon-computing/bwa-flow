#ifndef OPENCLENV_H
#define OPENCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <deque>
#include <glog/logging.h>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

#include "config.h"
#include "kflow/Queue.h"
#include "util.h"

#define OCL_CHECK(err, msg) { \
    if (err != CL_SUCCESS) { \
      throw std::runtime_error(msg); \
    } \
  }

#define OCL_CHECKRUN(func_call, msg) { \
    cl_int err = (func_call); \
    OCL_CHECK(err, msg); \
  }

struct cl_device_env {
  cl_context       context;
  cl_device_id     device_id;
  cl_command_queue cmd;
  cl_program       program;
  cl_device_env(): context(NULL), device_id(NULL), cmd(NULL), program(NULL) {;}
};

const cl_device_env NULL_DEVICE_ENV;

class OpenCLEnv 
: public boost::basic_lockable_adapter<boost::mutex> {
  public:
    OpenCLEnv(const char* bin_path, const char* kernel_name) {
      // start platform setting up
      int err;

      cl_platform_id* platform_id = new cl_platform_id[2];

      char cl_platform_vendor[1001];
      char cl_platform_name[1001];

      cl_uint num_platforms = 0;

      // Connect to first platform
      OCL_CHECKRUN(clGetPlatformIDs(2, platform_id, &num_platforms),
          "Failed to find an OpenCL platform!");

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

        if (strstr(cl_platform_name, "Xilinx") != NULL || 
            strstr(cl_platform_name, "Intel") != NULL ||
            strstr(cl_platform_name, "Altera") != NULL) {
          // found platform
          platform_idx = i;
          DLOG(INFO) << "Selecting platform " << cl_platform_name;
          break;
        }
      }

      // Load binary from disk
      unsigned char *kernelbinary;

      int n_i = load_file(bin_path, (char **) &kernelbinary);
      if (n_i < 0) {
        throw std::runtime_error(
            "failed to load kernel from xclbin");
      }
      size_t n_t = n_i;

      // Connect to a compute devices
      cl_uint num_devices = 0;
      int max_devices = FLAGS_max_fpga_thread;
      if (max_devices <= 0) {
        throw std::runtime_error("No FPGA thread configured,"
            " do not create platform");
      }

      cl_device_id* device_ids = (cl_device_id*)malloc(
          max_devices*sizeof(cl_device_id));

      err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR,
          max_devices, device_ids, &num_devices);

      if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group: " +
            std::to_string((long long)err));
      }
      DLOG(INFO) << "Found " << num_devices << " devices";

      for (int i = 0; i < num_devices; i++) {
        cl_int err = 0;
        int status = 0;

        cl_device_env env;

        env.device_id = device_ids[i];
        env.context = clCreateContext(0, 1, &device_ids[i], NULL, NULL, &err);
        OCL_CHECK(err, "failed to create context");

        env.cmd = clCreateCommandQueue(env.context, env.device_id, 0, &err);
        OCL_CHECK(err, "failed to create cmd_queue");

        env.program = clCreateProgramWithBinary(env.context, 1, &env.device_id, &n_t,
            (const unsigned char **) &kernelbinary, &status, &err);
        OCL_CHECK(err, "failed to create program from binary");

        OCL_CHECKRUN(clBuildProgram(env.program, 0, NULL, NULL, NULL, NULL), 
            "failed to build program executable");
        device_envs_.push_back(env);
      }

      free(device_ids);
    }

    ~OpenCLEnv() {
      // may not obtain all the devices back
      // need to pay attention to the class destroy
      // order
      while (!device_envs_.empty()) {
        cl_device_env env = device_envs_.front();

        clReleaseCommandQueue(env.cmd);
        clReleaseContext(env.context);
        clReleaseProgram(env.program);

        device_envs_.pop_front();
      }
    }

    // exclusive access to a device's context and queue
    cl_device_env getDevice() {
      boost::lock_guard<OpenCLEnv> guard(*this);
      if (device_envs_.empty()) {
        DLOG(ERROR) << "No more device available on this platform";
        return NULL_DEVICE_ENV;
      }
      // return next device and accumulate the index
      cl_device_env ret = device_envs_.front();
      device_envs_.pop_front();
      return ret;
    }

    void releaseDevice(cl_device_env env) {
      boost::lock_guard<OpenCLEnv> guard(*this);
      device_envs_.push_back(env);
      DLOG(INFO) << "release one opencl device";
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

  protected:
    std::deque<cl_device_env> device_envs_;
};

#endif

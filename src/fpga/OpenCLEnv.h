#ifndef OPENCLENV_H
#define OPENCLENV_H

#include <boost/atomic.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/lockable_adapter.hpp>
#include <boost/thread.hpp>
#include <deque>
#include <glog/logging.h>
#include <map>
#include <string>
#include <stdexcept>

#include <CL/opencl.h>

#include "config.h"
#include "kflow/Queue.h"
#include "util.h"
#include "allocation_wrapper.h"

#define OCL_CHECK(err, msg) { \
    if (err != CL_SUCCESS) { \
      DLOG(ERROR) << msg << " (error code: " << err << ")"; \
      throw std::runtime_error(msg); \
    } \
  }

#define OCL_CHECKRUN(func_call, msg) { \
    cl_int err = (func_call); \
    OCL_CHECK(err, msg); \
  }

struct cl_accx {
  int              accx_id;
  int              accx_group_id;
  cl_device_id     device_id;
  cl_context       context;
  cl_command_queue cmd;
  cl_program       program;
  cl_accx(): accx_id(-1), accx_group_id(-1), context(NULL), device_id(NULL), cmd(NULL), program(NULL) {;}
};

const cl_accx NULL_DEVICE_ENV;

class OpenCLEnv 
: public boost::basic_lockable_adapter<boost::mutex> 
{
  public:
    OpenCLEnv(std::vector<int> &group_sizes, std::vector<const char*> &bin_paths): device_envs_()
    {
      boost::lock_guard<OpenCLEnv> guard(*this);
      // start platform setting up
      cl_int err;

      cl_platform_id platform_id[2];

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

      // Split devices into groups 
      cl_uint num_devices = getMaxNumDevices();

      cl_device_id* device_ids = (cl_device_id*)malloc(
          num_devices*sizeof(cl_device_id));

      err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR,
          num_devices, device_ids, NULL);
      if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to create a device group: " +
            std::to_string((long long)err));
      }

      int max_devices = FLAGS_max_fpga_thread;

      if (FLAGS_max_fpga_thread < 0) {
        max_devices = num_devices;
      } else {
        max_devices = std::min(FLAGS_max_fpga_thread, int(num_devices));
      }
      DLOG(INFO) << "Found " << num_devices << " devices, "
                 << "using " << max_devices;

      int num_groups = std::min(group_sizes.size(), bin_paths.size());
      //int total_num = 0;
      //for (int i = 0; i < num_groups; i++) {
      //  if (total_num < 0 || group_sizes[i] < 0)
      //    total_num = -1;
      //  else
      //    total_num += group_sizes[i];
      //}
      //if (total_num < 0 || total_num > max_devices) {
      //  for (int i = 0; i < num_groups; i++) {
      //     group_sizes[i] = max_devices/num_groups;
      //     if (i < max_devices%num_groups) group_sizes[i] += 1;
      //  }
      //}

      int total_id = 0;
      for (int gid = 0; gid < num_groups; gid++) {
        if (group_sizes[gid] <= 0) continue;

        // Load binary from disk
        unsigned char *kernelbinary;

        int n_i = load_file(bin_paths[gid], (char **) &kernelbinary);
        if (n_i < 0) {
          throw std::runtime_error("failed to load kernel from xclbin - " + std::string(bin_paths[gid]));
        }
        size_t n_t = n_i;

        for (int i = 0; i < group_sizes[gid]; i++) {
          cl_int err = 0;
          cl_int status = 0;

          cl_accx env;

          env.accx_id = total_id;
          env.accx_group_id = gid;
          env.device_id = device_ids[total_id];
          DLOG(INFO) << "creating context for device #" << total_id;

          env.context = clCreateContext(NULL, 1, &env.device_id, NULL, NULL, &err);
          OCL_CHECK(err, "failed to create context");

          //env.cmd = clCreateCommandQueue(env.context, env.device_id, CL_QUEUE_PROFILING_ENABLE|CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
#ifdef XILINX_FPGA
          env.cmd = clCreateCommandQueue(env.context, env.device_id, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, &err);
#else
          env.cmd = clCreateCommandQueue(env.context, env.device_id, 0, &err);
#endif
          OCL_CHECK(err, "failed to create cmd_queue");

          env.program = clCreateProgramWithBinary(env.context, 1, &env.device_id, &n_t,
              (const unsigned char **) &kernelbinary, &status, &err);
          if (err == CL_INVALID_CONTEXT) DLOG(ERROR) << "ctx";
          if (err == CL_INVALID_VALUE) DLOG(ERROR) << "val";
          if (err == CL_INVALID_DEVICE) DLOG(ERROR) << "dvc";
          if (err == CL_INVALID_BINARY) DLOG(ERROR) << "bin";
          if (err == CL_OUT_OF_HOST_MEMORY) DLOG(ERROR) << "mem";
          OCL_CHECK(err, "failed to create program from binary");

          OCL_CHECKRUN(clBuildProgram(env.program, 0, NULL, NULL, NULL, NULL), 
              "failed to build program executable");

          device_envs_.push_back(env);
          DLOG(INFO) << "setup opencl env for one device";
          total_id++;
        }
      }

      FLAGS_max_fpga_thread = max_devices;
      free(device_ids);
    }

    ~OpenCLEnv() {
      // may not obtain all the devices back
      // need to pay attention to the class destroy
      // order
      while (!device_envs_.empty()) {
        cl_accx env = device_envs_.front();

        clReleaseCommandQueue(env.cmd);
        clReleaseContext(env.context);
        clReleaseProgram(env.program);

        device_envs_.pop_front();
      }
    }

    int getMaxNumDevices() {
      cl_int err;

      cl_platform_id platform_id[2];

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

      cl_uint num_devices = 0;

      err = clGetDeviceIDs(platform_id[platform_idx], CL_DEVICE_TYPE_ACCELERATOR,
          0, NULL, &num_devices);
      if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to request the devices number: " +
            std::to_string((long long)err));
      }

#ifdef DEPLOY_aws
      //num_devices = num_devices/2;
#endif

      return num_devices;
    }

/*
    // exclusive access to a device's context and queue
    cl_accx getDevice() {
      boost::lock_guard<OpenCLEnv> guard(*this);
      uint32_t tid = getTid();
      if (device_registry_.count(tid)) {
        // already allocated a device to this thread
        DLOG(INFO) << "return allocated device for " << tid;
        device_ref_counter_[tid]++;
        return device_registry_[tid];
      }
      else if (device_envs_.empty()) {
        DLOG(ERROR) << "No more device available on this platform";
        return NULL_DEVICE_ENV;
      }
      else {
        // return next device and accumulate the index
        cl_accx ret = device_envs_.front();
        device_envs_.pop_front();
        device_registry_[tid] = ret;
        device_ref_counter_[tid] = 0;
        DLOG(INFO) << "allocate one opencl device for " << tid;
        return ret;
      }
    }

    void releaseDevice(cl_accx env) {
      boost::lock_guard<OpenCLEnv> guard(*this);
      uint32_t tid = getTid();
      device_ref_counter_[tid]--;
      if (device_ref_counter_[tid] == 0) {
        device_registry_.erase(tid);
        device_ref_counter_.erase(tid);
        device_envs_.push_back(env);
      }
      DLOG(INFO) << "release one opencl device";
    }
*/

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

/*
    // get current thread id
    // using the same code from googlelog/src/utilities.cc
    // without OS checking
    uint32_t getTid() {
      static bool lacks_gettid = false;

      if (!lacks_gettid) {
        pid_t tid = syscall(__NR_gettid);
        if (tid != -1) {
          return (uint32_t)tid;
        }
        // Technically, this variable has to be volatile, but there is a small
        // performance penalty in accessing volatile variables and there should
        // not be any serious adverse effect if a thread does not immediately see
        // the value change to "true".
        lacks_gettid = true;
      }

      // If gettid() could not be used, we use one of the following.
      return (uint32_t)getpid();
    }
*/

  protected:
    std::map<uint32_t, cl_accx> device_registry_;
    std::map<uint32_t, int> device_ref_counter_;
    std::deque<cl_accx> device_envs_;
};

#endif

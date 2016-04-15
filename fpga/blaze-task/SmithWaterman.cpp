#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>

#include "blaze/Task.h" 
#include "OpenCLEnv.h" 

#define FPGA_RET_PARAM_NUM 5

using namespace blaze;

class SmithWaterman : public Task {
public:

  // extends the base class constructor
  // to indicate how many input blocks
  // are required
  SmithWaterman(): Task(2) {; }

  virtual void compute() {

    struct	timeval t1, t2, tr;

    try {
      // dynamically cast the TaskEnv to OpenCLEnv
      OpenCLEnv* ocl_env = (OpenCLEnv*)getEnv();

      // get input data length
      int data_size = getInputLength(0);
      int task_num  = *((int*)getInput(1));

      // check input size
      if (data_size < 1 || task_num < 1 ) 
      {
        throw std::runtime_error("Invalid input data dimensions");
      }

      // get OpenCL context
      cl_context       context = ocl_env->getContext();
      cl_kernel        kernel  = ocl_env->getKernel();
      cl_command_queue command = ocl_env->getCmdQueue();

      int err;
      cl_event event;

      // get the pointer to input/output data
      cl_mem ocl_input  = *((cl_mem*)getInput(0));
      cl_mem ocl_output = *((cl_mem*)getOutput(
            0, FPGA_RET_PARAM_NUM*task_num, 1, sizeof(int)));

      if (!ocl_input || !ocl_output) {
        throw std::runtime_error("Buffer are not allocated");
      }

      // Set the arguments to our compute kernel
      err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &ocl_input);
      err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &ocl_output);
      err |= clSetKernelArg(kernel, 2, sizeof(int), &task_num);
      if (err != CL_SUCCESS) {
        throw std::runtime_error("Failed to set gradients!");
      }

      // Execute the kernel over the entire range of our 1d input data set
      // using the maximum number of work group items for this device
      ocl_env->lock();

      err = clEnqueueTask(command, kernel, 0, NULL, &event);
      if (err) {
        throw("Failed to execute kernel!");
      }

      ocl_env->unlock();
      clWaitForEvents(1, &event);
    }
    catch (std::runtime_error &e) {
      throw e;
    }
  }
};

extern "C" Task* create() {
  return new SmithWaterman();
}

extern "C" void destroy(Task* p) {
  delete p;
}

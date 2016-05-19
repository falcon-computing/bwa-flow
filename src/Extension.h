#ifndef EXTENSION_H
#define EXTENSION_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "FPGAAgent.h"
#include "SWTask.h"

extern OpenCLEnv* opencl_env;

extern  ExtParam** task_batch[5];
extern  bool pend_flag[5];
extern  int pend_depth; 
extern  boost::mutex driver_mutex;
extern  boost::condition_variable driver_cond;

void extendOnFPGAPackInput(
    FPGAAgent* agent,
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt);

void extendOnFPGAProcessOutput(
    FPGAAgent* agent,
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt);

void fpga_driver(
    FPGAAgent* agent);

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt);

#endif

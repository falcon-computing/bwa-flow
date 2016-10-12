#ifndef EXTENSION_H
#define EXTENSION_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "FPGAAgent.h"
#include "SWTask.h"

extern OpenCLEnv* opencl_env;

void extendOnFPGA(
    FPGAAgent* agent,
    char* &kernel_buffer,
    int data_size,
    int stage_cnt
    );

void FPGAPostProcess(
    FPGAAgent* agent,
    short* kernel_output,
    int task_num,
    mem_alnreg_t** &region_batch,
    mem_chain_t** &chain_batch,
    int stage_cnt
    );
#endif

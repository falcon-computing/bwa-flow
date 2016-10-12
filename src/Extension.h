#ifndef EXTENSION_H
#define EXTENSION_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "FPGAAgent.h"
#include "SWTask.h"

extern OpenCLEnv* opencl_env;

void buffer_out_process(short* kernel_output, int task_num, mem_alnreg_t** &region_batch);

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
    int stage_cnt
    );

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

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt);

#endif

#ifndef EXTENSION_H
#define EXTENSION_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "FPGAAgent.h"
#include "SWTask.h"

extern OpenCLEnv* opencl_env;

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

char* extendOnFPGAonlyPack(
     ExtParam** &tasks,
     int batch_num,
     mem_opt_t *opt,
     int *fpga_data_size);

void extendOnFPGAonlyOutput (
     ExtParam** &tasks,
     int batch_num,
     short* &output_ptr,
     mem_opt_t* opt);

#endif

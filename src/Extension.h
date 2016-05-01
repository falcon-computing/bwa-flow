#ifndef EXTENSION_H
#define EXTENSION_H

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "FPGAAgent.h"
#include "SWTask.h"

extern FPGAAgent* agent;

void packData(int stage_cnt,
      ExtParam** &tasks,
      int batch_num,
      mem_opt_t *opt);

void SwFPGA(
    int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt);

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt);

#endif

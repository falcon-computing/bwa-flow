#ifndef EXTENSION_H
#define EXTENSION_H

#include "SWTask.h"

#define MAX_BAND_TRY  2

blaze::Task_ptr packData(
      ExtParam** &tasks,
      int batch_num,
      mem_opt_t *opt);

void SwFPGA(
    ExtParam** &tasks,
    blaze::Task_ptr fpga_task,
    int batch_num,
    mem_opt_t *opt);

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt);

#endif

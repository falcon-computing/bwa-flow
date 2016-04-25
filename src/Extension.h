#ifndef EXTENSION_H
#define EXTENSION_H

#include "SWTask.h"

#define MAX_BAND_TRY  2

void packData(int stage_cnt,
      ExtParam** &tasks,
      int batch_num,
      mem_opt_t *opt);

void SwFPGA(int stage_cnt,
    ExtParam** &tasks,
    int batch_num,
    mem_opt_t *opt);

void extendOnCPU(
    ExtParam** tasks,
    int numoftask,
    mem_opt_t *opt);

#endif

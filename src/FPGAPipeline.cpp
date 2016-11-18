#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <vector>
#include <queue>
#include <list>
#include <boost/asio.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>

#include "bwa_wrapper.h"
#include "config.h"
#include "FPGAAgent.h"
#include "util.h"
#include "Pipeline.h"

#define MAX_BAND_TRY  2
//#define SMITHWATERMAN_SIM

#ifdef SMITHWATERMAN_SIM
// hw data structures
extern "C"{
void sw_top (int *a, int *output_a, int __inc);
}
#endif

void ChainsToRegionsFPGA::extendOnFPGA(
    FPGAAgent* agent,
    char* &kernel_buffer,
    int data_size_a,
    int data_size_b,
    int stage_cnt
    )
{
  agent->writeInput((int*)kernel_buffer, data_size_a, stage_cnt, 0);
  agent->writeInput((int*)(&kernel_buffer[data_size_a]), data_size_b, stage_cnt, 1);
  agent->start(data_size_a/4, data_size_b/4, stage_cnt);
}

void ChainsToRegionsFPGA::FPGAPostProcess(
    FPGAAgent* agent,
    short* kernel_output,
    int task_num_a,
    int task_num_b,
    mem_alnreg_t** &region_batch,
    mem_chain_t** &chain_batch,
    int stage_cnt
    )
{
  uint64_t start_ts = getUs();
  int seedcov = 0;
  agent->wait(stage_cnt);
  VLOG(3) << "Wait for FPGA takes " 
    << getUs() - start_ts << " us";

  start_ts = getUs();
  agent->readOutput(kernel_output, FPGA_RET_PARAM_NUM*task_num_a*4, stage_cnt, 0);
  agent->readOutput(&kernel_output[FPGA_RET_PARAM_NUM*task_num_a*2],
      FPGA_RET_PARAM_NUM*task_num_b*4, stage_cnt, 1);
  for (int i = 0; i < task_num_a + task_num_b; i++) {
    int seed_idx = ((int)(kernel_output[1+FPGA_RET_PARAM_NUM*2*i])<<16) |
                    kernel_output[0+FPGA_RET_PARAM_NUM*2*i];
    mem_alnreg_t *newreg = region_batch[seed_idx];

    newreg->qb = kernel_output[2+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->qe += kernel_output[3+FPGA_RET_PARAM_NUM*2*i];
    newreg->rb += kernel_output[4+FPGA_RET_PARAM_NUM*2*i];
    newreg->re += kernel_output[5+FPGA_RET_PARAM_NUM*2*i];
    newreg->score = kernel_output[6+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->truesc = kernel_output[7+FPGA_RET_PARAM_NUM*2*i]; 
    newreg->w = kernel_output[8+FPGA_RET_PARAM_NUM*2*i];

    // compute the seedcov here
    mem_chain_t *chain = chain_batch[seed_idx];
    seedcov = 0;
    for (int k=0; k < chain->n; k++) {
      const mem_seed_t *seed = &chain->seeds[k];
      if (seed->qbeg >= newreg->qb && 
          seed->qbeg + seed->len <= newreg->qe && 
          seed->rbeg >= newreg->rb && 
          seed->rbeg + seed->len <= newreg->re){
        seedcov += seed->len; 
      }
    }
    newreg->seedcov = seedcov;
  }
  VLOG(3) << "Process output takes " 
    << getUs() - start_ts << " us";
}

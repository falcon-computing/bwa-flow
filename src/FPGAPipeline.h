#ifndef BWA_FLOW_FPGAPIPELINE_H
#define BWA_FLOW_FPGAPIPELINE_H

#include "bwa_wrapper.h"
#include "Pipeline.h"

#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"

class ChainsPipeFPGA
: public kestrelFlow::MapStage<
      ChainsRecord, ChainsRecord, COMPUTE_DEPTH, 8>
{
 public:
  ChainsPipeFPGA(int n=1): kestrelFlow::MapStage<
      ChainsRecord, ChainsRecord, COMPUTE_DEPTH, 8>(n, false) {;}

  ChainsRecord compute(ChainsRecord const & record) {
    ChainsRecord output = record;
    return output; 
  }
};

class ChainsToRegionsFPGA
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, 8, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsFPGA(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, 8, COMPUTE_DEPTH>(n, false) {;}

  void compute(int wid);
  void extendOnFPGA(
      FPGAAgent* agent,
      char* &kernel_buffer,
      int data_size_a,
      int data_size_b,
      int stage_cnt
      );
  void FPGAPostProcess(
      FPGAAgent* agent,
      short* kernel_output,
      int task_num_a,
      int task_num_b,
      mem_alnreg_t** &region_batch,
      mem_chain_t** &chain_batch,
      int stage_cnt
      );
};

#endif


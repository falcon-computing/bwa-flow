#ifndef BWA_FLOW_FPGAPIPELINE_H
#define BWA_FLOW_FPGAPIPELINE_H

#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "SWTask.h"

#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"


class ChainsToRegionsFPGA
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsFPGA(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  void compute(int wid);
  void processOutput(SWTask* task, uint64_t &deq_ts, uint64_t &read_ts, uint64_t &post_ts);
};

class SeqsToChainsFPGA
: public kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  SeqsToChainsFPGA(int n=1) : kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  //ChainsRecord compute(SeqsRecord const & record);
  void compute(int wid);
};

#endif


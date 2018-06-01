#ifndef BWA_FLOW_FPGAPIPELINE_H
#define BWA_FLOW_FPGAPIPELINE_H

#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "SWTask.h"

#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"

class ChainsPipeFPGA
: public kestrelFlow::MapStage<
      ChainsRecord, ChainsRecord, COMPUTE_DEPTH, 8>
{
 public:
  ChainsPipeFPGA(int n=1): kestrelFlow::MapStage<
      ChainsRecord, ChainsRecord, COMPUTE_DEPTH, 8>(n, false) {;}

  ChainsRecord compute(ChainsRecord const & record) ;
};

class ChainsToRegionsFPGA
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, 8, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsFPGA(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, 8, COMPUTE_DEPTH>(n, false) {;}

  void compute(int wid);
  void processOutput(SWTask* task);
};

class SeqsToChainsFPGA
: public kestrelFlow::MapStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  SeqsToChainsFPGA(int n=1) : kestrelFlow::MapStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  ChainsRecord compute(SeqsRecord const & record);
};

#endif


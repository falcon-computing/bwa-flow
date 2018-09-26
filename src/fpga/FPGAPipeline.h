#ifndef BWA_FLOW_FPGAPIPELINE_H
#define BWA_FLOW_FPGAPIPELINE_H

#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "SWTask.h"
#include "SMemTask.h"

#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "boost/atomic.hpp"
#include "boost/thread/mutex.hpp"

class ChainsToRegionsFPGA
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsFPGA(int n=1, ChainsToRegions *stage=NULL) : kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false),
      n_active_(n),
      cpu_stage_(stage)
  {;}

  void compute(int wid);
  void swtask_func(int wid, void *param_list[3]);

 private:
  boost::atomic<int> n_active_;
  ChainsToRegions * cpu_stage_;
};

class SeqsToChainsFPGA
: public kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  SeqsToChainsFPGA(int n=1, SeqsToChains *stage=NULL) : kestrelFlow::MapPartitionStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH>(n, false),
      n_active_(n),
      cpu_stage_(stage)
  {;}

  void compute(int wid);

 private:
  boost::atomic<int> n_active_;
  SeqsToChains * cpu_stage_;
};

#endif


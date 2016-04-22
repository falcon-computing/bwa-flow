#ifndef BWA_FLOW_PIPELINE_H
#define BWA_FLOW_PIPELINE_H

#include <list>
#include <unordered_map>

#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "bwa_wrapper.h"
#include "SWRead.h"

// Common data structures
struct SeqsRecord {
  int start_idx;
  int batch_num;
  bseq1_t* seqs;
};

struct ChainsRecord {
  int start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
  std::list<SWRead*>* read_batch;
};

struct RegionsRecord {
  int start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_alnreg_v* alnreg;
};

class SeqsProducer : public kestrelFlow::SourceStage<SeqsRecord> {
 public:
  SeqsProducer(): kestrelFlow::SourceStage<SeqsRecord>() {;}
  void compute();
};

class SeqsToChains 
: public kestrelFlow::MapStage<SeqsRecord, ChainsRecord> {
 public:
  SeqsToChains(int n=1): 
      kestrelFlow::MapStage<SeqsRecord, ChainsRecord>(n),
      aux(NULL) {;}

  ChainsRecord compute(SeqsRecord const & record);
 private:
  ktp_aux_t* aux;
};

// TODO: this in the future may not be a map stage anymore
class ChainsToRegions
: public kestrelFlow::MapPartitionStage<ChainsRecord, RegionsRecord> {
 public:
  ChainsToRegions(int n=1): 
      kestrelFlow::MapPartitionStage<ChainsRecord, RegionsRecord>(n) {;}

  void compute();
 private:
  bool addBatch();

  // Batch of SWTasks
  ExtParam**         task_batch_;

  // Batch of SWReads
  std::list<SWRead*> read_batch_;

  // Table to keep track of each record
  std::unordered_map<uint64_t, int> tasks_remain_;
  std::unordered_map<uint64_t, ChainsRecord> input_buf_;
  std::unordered_map<uint64_t, RegionsRecord> output_buf_;
};

class RegionsToSam : public kestrelFlow::MapStage<RegionsRecord, SeqsRecord> {
 public:
  RegionsToSam(int n=1): 
      kestrelFlow::MapStage<RegionsRecord, SeqsRecord>(n),
      aux(NULL) {;}

  SeqsRecord compute(RegionsRecord const & record);
 private:
  ktp_aux_t* aux;
};

class PrintSam : public kestrelFlow::SinkStage<SeqsRecord> {
 public:
  PrintSam(): kestrelFlow::SinkStage<SeqsRecord>() {;}
  void compute();
};

#endif

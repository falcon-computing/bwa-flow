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
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
};

struct ChainsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
  std::vector<int>* chains_idxes;
  std::list<SWRead*>* read_batch;
};

struct RegionsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
  std::vector<int>* chains_idxes;
};

class SeqsProducer : public kestrelFlow::SourceStage<SeqsRecord, 4> {
 public:
  SeqsProducer(): kestrelFlow::SourceStage<SeqsRecord, 4>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
  SeqsRecord deserialize(const char* data, size_t length);
};

class SeqsToChains 
: public kestrelFlow::MapStage<SeqsRecord, ChainsRecord, 4, 16> {
 public:
  SeqsToChains(int n=1): 
      kestrelFlow::MapStage<SeqsRecord, ChainsRecord, 4, 16>(n),
      aux(NULL) {;}

  ChainsRecord compute(SeqsRecord const & record);
 private:
  ktp_aux_t* aux;
};

class ChainsToRegions
: public kestrelFlow::MapPartitionStage<ChainsRecord, RegionsRecord, 16, 16>
{
 public:
  ChainsToRegions(int n=1): 
      kestrelFlow::MapPartitionStage<ChainsRecord, RegionsRecord, 16, 16>(n) {;}

  void compute(int wid);
 private:
  inline bool addBatch(
      std::list<SWRead*> &read_batch,
      std::unordered_map<uint64_t, int> &tasks_remain,
      std::unordered_map<uint64_t, ChainsRecord> &input_buf,
      std::unordered_map<uint64_t, RegionsRecord> &output_buf);

};

class RegionsToSam
: public kestrelFlow::MapStage<RegionsRecord, SeqsRecord, 16, 4> 
{
 public:
  RegionsToSam(int n=1): 
      kestrelFlow::MapStage<RegionsRecord, SeqsRecord, 16, 4>(n),
      aux(NULL) {;}

  SeqsRecord compute(RegionsRecord const & record);
 private:
  ktp_aux_t* aux;
};

class PrintSam : public kestrelFlow::SinkStage<SeqsRecord, 4> {
 public:
  PrintSam(): kestrelFlow::SinkStage<SeqsRecord, 4>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
  SeqsRecord deserialize(const char* data, size_t length);
};

#endif

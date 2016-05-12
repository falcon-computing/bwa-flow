#ifndef BWA_FLOW_PIPELINE_H
#define BWA_FLOW_PIPELINE_H

#include <boost/thread/mutex.hpp>
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

#define INPUT_DEPTH   4
#define OUTPUT_DEPTH  8

#define MASTER_RANK   0

// Tags for messages between master and child processes
#define SEQ_DP_QUERY  0
#define SEQ_DP_LENGTH 1
#define SEQ_DP_DATA   2
#define SAM_RV_QUERY  3
#define SAM_RV_LENGTH 4
#define SAM_RV_DATA   5

// mutex for serializing MPI calls
extern boost::mutex mpi_mutex;

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

#ifdef SCALE_OUT
class SeqsDispatch : public kestrelFlow::SinkStage<SeqsRecord, 4> {
 public:
  SeqsDispatch(): kestrelFlow::SinkStage<SeqsRecord, 4>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
};

class SeqsReceive : public kestrelFlow::SourceStage<SeqsRecord, 4> {
 public:
  SeqsReceive(): kestrelFlow::SourceStage<SeqsRecord, 4>() {;}
  void compute();
  SeqsRecord deserialize(const char* data, size_t length);
};

class SamsSend : public kestrelFlow::SinkStage<SeqsRecord, 8> {
 public:
  SamsSend(): kestrelFlow::SinkStage<SeqsRecord, 8>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
};

class SamsReceive : public kestrelFlow::SourceStage<SeqsRecord, 8> {
 public:
  SamsReceive(): kestrelFlow::SourceStage<SeqsRecord, 8>() {;}
  SeqsRecord deserialize(const char* data, size_t length);
  void compute();
};
#endif

class SeqsRead : public kestrelFlow::SourceStage<SeqsRecord, 4> {
 public:
  SeqsRead(): kestrelFlow::SourceStage<SeqsRecord, 4>() {;}
  void compute();
};

// One stage for the entire BWA-MEM computation
class SeqsToSams
: public kestrelFlow::MapStage<SeqsRecord, SeqsRecord, 4, 8> {
 public:
  SeqsToSams(int n=1): 
      kestrelFlow::MapStage<SeqsRecord, SeqsRecord, 4, 8>(n)
  {;}

  SeqsRecord compute(SeqsRecord const & record);
};

class SeqsToChains 
: public kestrelFlow::MapStage<SeqsRecord, ChainsRecord, 4, 16> {
 public:
  SeqsToChains(int n=1): 
      kestrelFlow::MapStage<SeqsRecord, ChainsRecord, 4, 16>(n)
  {;}

  ChainsRecord compute(SeqsRecord const & record);
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
: public kestrelFlow::MapStage<RegionsRecord, SeqsRecord, 16, 8> 
{
 public:
  RegionsToSam(int n=1): 
      kestrelFlow::MapStage<RegionsRecord, SeqsRecord, 16, 8>(n)
  {;}

  SeqsRecord compute(RegionsRecord const & record);
};

class SamsPrint : public kestrelFlow::SinkStage<SeqsRecord, 8> {
 public:
  SamsPrint(): kestrelFlow::SinkStage<SeqsRecord, 8>() {;}
  void compute();
};

#endif

#ifndef BWA_FLOW_PIPELINE_H
#define BWA_FLOW_PIPELINE_H

#include <boost/thread/mutex.hpp>
#include <glog/logging.h>
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
#include "FPGAAgent.h"

#define INPUT_DEPTH   4
#define OUTPUT_DEPTH  16
#define COMPUTE_DEPTH 16

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
};

struct RegionsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
};

#ifdef SCALE_OUT
class SeqsDispatch : public kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsDispatch(): kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
};

class SeqsReceive : public kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsReceive(): kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute();
  SeqsRecord deserialize(const char* data, size_t length);
};

class SamsSend : public kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsSend(): kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>() {;}
  void compute();
  std::string serialize(SeqsRecord* data);
};

class SamsReceive : public kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsReceive(): kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH>() {;}
  SeqsRecord deserialize(const char* data, size_t length);
  void compute();
};
#endif

class SeqsRead : public kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsRead(): kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute();
};

// One stage for the entire BWA-MEM computation
class SeqsToSams
: public kestrelFlow::MapStage<
    SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH> {
 public:
  SeqsToSams(int n=1): 
      kestrelFlow::MapStage<
          SeqsRecord, SeqsRecord, INPUT_DEPTH, OUTPUT_DEPTH>(n)
  {;}

  SeqsRecord compute(SeqsRecord const & record);
};

class SeqsToChains 
: public kestrelFlow::MapStage<
      SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH> {
 public:
  SeqsToChains(int n=1): 
      kestrelFlow::MapStage<
          SeqsRecord, ChainsRecord, INPUT_DEPTH, COMPUTE_DEPTH>(n)
  {;}
  ChainsRecord compute(SeqsRecord const & record);
};

class ChainsToRegions
: public kestrelFlow::MapStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegions(int n=1): kestrelFlow::MapStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;}

  RegionsRecord compute(ChainsRecord const & record);
};

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

class RegionsToSam
: public kestrelFlow::MapStage<
      RegionsRecord, SeqsRecord, COMPUTE_DEPTH, OUTPUT_DEPTH> 
{
 public:
  RegionsToSam(int n=1): kestrelFlow::MapStage<
      RegionsRecord, SeqsRecord, COMPUTE_DEPTH, OUTPUT_DEPTH>(n)
  {;}

  SeqsRecord compute(RegionsRecord const & record);
};

class SamsPrint
: public kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsPrint():
    kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>(),
    file_id_(0), fout_(NULL)
  {;}
  void compute();
 private:
  void sortAndWriteBamBatch(bam1_t** buf, int n_elements, std::string out_dir);

  int file_id_;
#ifdef USE_HTSLIB
  samFile* fout_;
#else
  FILE* fout_
#endif
};

#endif

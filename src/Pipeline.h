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
#include "SWRead.h"

#define INPUT_DEPTH   4
#define OUTPUT_DEPTH  16
#define COMPUTE_DEPTH 64

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
  std::list<SWRead*>* read_batch;
  std::vector<int>* chains_idxes;
  mem_chainref_t** chainrefs;
  char* fpga_data_buf;
  int fpga_data_size;
  int fpga_task_num;
  ExtParam** fpga_task_batch;
};

struct RegionsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
  std::vector<int>* chains_idxes;
  mem_chainref_t** chainrefs;
  short* fpga_output_buf;
  int fpga_task_num;
  ExtParam** fpga_task_batch;
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
 private:
  static inline void prepareChainRef(
      const ktp_aux_t* aux,
      const bseq1_t* seq,
      const mem_chain_v* chain,
      mem_chainref_t* &ref);
};

class ChainsToRegions
: public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegions(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n) {;}

  void compute(int wid);
 private:
  inline bool addBatch(
      std::list<SWRead*> &read_batch,
      std::unordered_map<uint64_t, int> &tasks_remain,
      std::unordered_map<uint64_t, ChainsRecord> &input_buf,
      std::unordered_map<uint64_t, RegionsRecord> &output_buf);

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
  
  inline int testExtension(
    mem_opt_t *opt,
    mem_seed_t& seed,
    mem_alnreg_v& alnregv, 
    int l_query); 
  inline int checkOverlap(
    int startidx,
    mem_seed_t& seed,
    mem_chain_t& chain,
    uint64_t *srt);
  inline void regionFilter(mem_alnreg_v alnreg,mem_alnreg_v &alnreg_short,
                mem_chainref_t* &chainref,mem_chain_v* chain,
                std::vector<int> &chain_idxes,
                bseq1_t* seq);
  
};

class SamsPrint
: public kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsPrint(): kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>() {;}
  void compute();
};

#endif

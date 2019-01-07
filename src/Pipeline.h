#ifndef BWA_FLOW_PIPELINE_H
#define BWA_FLOW_PIPELINE_H

#include <boost/thread/mutex.hpp>
#include <glog/logging.h>
#include <list>
#include <unordered_map>

#include "BamFileBuffer.h"
#include "bwa/bwamem.h"
#include "bwa/bntseq.h"
#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "bwa_wrapper.h"

#define INPUT_DEPTH   64
#define OUTPUT_DEPTH  64
#define COMPUTE_DEPTH 64

void sort_bams(int size, bam1_t** buffer);

// Common data structures
struct kseq_buf_t {
  kseq_new_t* ks;
  int size;
};

struct KseqsRecord {
  uint64_t start_idx;
  int batch_num;
  kseq_buf_t ks_buffer;
  const char* name = "KseqsRecord";
};

struct SeqsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  const char* name = "SeqsRecord";
};

struct ChainsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  bwtintv_t** bwtintvs;
  size_t* bwtintv_nums;
  mem_alnreg_v* alnreg;
  mem_chainref_t** chain_ref;
  const char* name = "ChainsRecord";
  int tag;
};

struct RegionsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_chain_v* chains;
  mem_alnreg_v* alnreg;
  const char* name = "RegionsRecord";
};

#ifdef USE_HTSLIB
struct BamsRecord {
  int bam_buffer_order;
  bam1_t** bam_buffer;
  int bam_buffer_idx;
  std::vector<SeqsRecord>* records_list;
  const char* name = "BamsRecord";
};
#endif

//Data structure for sort-merge pipe
struct BamRecord {
  int id;
  int size;
  bam1_t ** bams;

  BamFileBuffer* fbuf;
};

template<typename Record>
inline void freeRecord(Record &record) {
  DLOG(INFO) << "Undefined record type: " << record.name;
}

template<>
inline void freeRecord(SeqsRecord &record) {
  freeSeqs(record.seqs, record.batch_num);
}

template<>
inline void freeRecord(ChainsRecord &record) {
  freeSeqs(record.seqs, record.batch_num);
  freeChains(record.chains, record.batch_num);
}

template<>
inline void freeRecord(RegionsRecord &record) {
  freeSeqs(record.seqs, record.batch_num);
  freeAligns(record.alnreg, record.batch_num);
}

class SeqsRead : public kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsRead(): kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute();
};

class KseqsRead : public kestrelFlow::SourceStage<KseqsRecord, INPUT_DEPTH> {
 public:
  KseqsRead(): kestrelFlow::SourceStage<KseqsRecord, INPUT_DEPTH>() {;}
  void compute();
};

class KseqsToBseqs : public kestrelFlow::MapStage<KseqsRecord, SeqsRecord, INPUT_DEPTH, INPUT_DEPTH> {
 public:
   KseqsToBseqs(int n=1):
     kestrelFlow::MapStage<KseqsRecord, SeqsRecord, INPUT_DEPTH, INPUT_DEPTH>(n) {;}
   SeqsRecord compute(KseqsRecord const &input);
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

class ChainsPipe
: public kestrelFlow::MapStage<
      ChainsRecord, ChainsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH> {
 public:
  ChainsPipe(int n=1): 
      kestrelFlow::MapStage<
          ChainsRecord, ChainsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n)
  {;}
  ChainsRecord compute(ChainsRecord const & record);
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


#ifdef USE_HTSLIB
class SamsReorder
: public kestrelFlow::MapPartitionStage<SeqsRecord, BamsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
  public:
    SamsReorder():kestrelFlow::MapPartitionStage
                  <SeqsRecord, BamsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(1, false){;}
    void compute(int wid);
};

class SamsSort
: public kestrelFlow::MapStage<BamsRecord, BamsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
  public:
    SamsSort(int n=1):kestrelFlow::MapStage
                  <BamsRecord, BamsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n){;}
    BamsRecord compute(BamsRecord const & record_bundle);
};
#else
class SamsReorder
: public kestrelFlow::MapPartitionStage<SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
  public:
    SamsReorder():kestrelFlow::MapPartitionStage
                  <SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(1, false){;}
    void compute(int wid);
};

class SamsSort
: public kestrelFlow::MapPartitionStage<SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
  public:
    SamsSort(int n=1):kestrelFlow::MapPartitionStage
                  <SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(1, false)
  {;}
    void compute(int wid);
};
#endif


class WriteOutput
#ifdef USE_HTSLIB
:public kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    WriteOutput(int n=1):
      kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, 0>(n, false) {;}
    int compute(BamsRecord const &input);
#else
:public kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, 0> {
  public:
    WriteOutput(int n=1):
      kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, 0>(n, false) {;}
    int compute(SeqsRecord const &input);
#endif
};

#endif

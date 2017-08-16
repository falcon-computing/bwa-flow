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

#define INPUT_DEPTH   8
#define OUTPUT_DEPTH  16
#define COMPUTE_DEPTH 16

// Common data structures
struct KseqsRecord {
  uint64_t start_idx;
  int batch_num;
  kseq_new_t* ks_buffer;
};

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
  mem_chainref_t** chain_ref;
};

struct RegionsRecord {
  uint64_t start_idx;
  int batch_num;
  bseq1_t* seqs;
  mem_alnreg_v* alnreg;
};

#ifdef USE_HTSLIB
struct BamsRecord {
  int bam_buffer_order;
  bam1_t** bam_buffer;
  int bam_buffer_idx;
};
#endif

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
                  <SeqsRecord, BamsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(1, false)
  {;}
    void compute(int wid);
};
#else
class SamsReorder
: public kestrelFlow::MapPartitionStage<SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
  public:
    SamsReorder():kestrelFlow::MapPartitionStage
                  <SeqsRecord, SeqsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(1, false)
  {;}
    void compute(int wid);
};
#endif


class WriteOutput
#ifdef USE_HTSLIB
:public kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, OUTPUT_DEPTH> {
  public:
    WriteOutput(int n=1):
      kestrelFlow::MapStage<BamsRecord, int, COMPUTE_DEPTH, OUTPUT_DEPTH>(n) {;}
    int compute(BamsRecord const &input);
#else
:public kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, OUTPUT_DEPTH> {
  public:
    WriteOutput(int n=1):
      kestrelFlow::MapStage<SeqsRecord, int, COMPUTE_DEPTH, OUTPUT_DEPTH>(n) {;}
    int compute(SeqsRecord const &input);
#endif
};

#endif

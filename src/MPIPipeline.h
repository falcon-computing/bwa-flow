#ifndef BWA_FLOW_MPIPIPELINE_H
#define BWA_FLOW_MPIPIPELINE_H

#define MASTER_RANK   0

// Tags for messages between master and child processes
#define SEQ_DP_QUERY  0
#define SEQ_DP_LENGTH 1
#define SEQ_DP_DATA   2
#define SAM_RV_QUERY  3
#define SAM_RV_LENGTH 4
#define SAM_RV_DATA   5

#include "Pipeline.h"
#include "MPIChannel.h"

// mutex for serializing MPI calls
extern boost::mutex mpi_mutex;

// MPI helper functions
namespace bwa_mpi {

// basic MPI APIs with mutex lock for multi-thread safety
bool query(MPI::Request &req);
void send(const void* buf, int count, const MPI::Datatype& datatype, 
    int dest, int tag);
void recv(void* buf, int count, const MPI::Datatype& datatype, 
    int source, int tag);

} // namepsace bwa_mpi

class SeqsDispatch : public kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsDispatch(MPILink* link): 
    ch_(link),
    kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH>() 
  {;}

  void compute(int wid = 0);

 private:
  SourceChannel ch_;
};

class SeqsReceive : public kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsReceive(MPILink* link): 
    ch_(link),
    kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH>() 
  {;}

  void compute();

 private:
  SourceChannel ch_;
};

/* 
 * Offload ChainsRecord to remote worker(s) and receive RegionsRecord.
 */
class ChainsToRegionsOffload : public kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>
{
 public:
  ChainsToRegionsOffload(int n=1): kestrelFlow::MapPartitionStage<
      ChainsRecord, RegionsRecord, COMPUTE_DEPTH, COMPUTE_DEPTH>(n, false) {;}

  void compute(int wid);
};

/*
 * Receive from remote process(es), and deserialize to ChainsRecord
 */
class ChainsReceive : public kestrelFlow::SourceStage<ChainsRecord, INPUT_DEPTH> {
 public:
  ChainsReceive(): kestrelFlow::SourceStage<ChainsRecord, INPUT_DEPTH>() {;}
  void compute();
};

/*
 * Serialize RegionsRecord and send to remote 
 */
class RegionsSend : public kestrelFlow::SinkStage<RegionsRecord, INPUT_DEPTH> {
 public:
  RegionsSend(): kestrelFlow::SinkStage<RegionsRecord, INPUT_DEPTH>() {;}
  void compute();
};

class SamsSend : public kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsSend(MPILink* link): 
    ch_(link),
    kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>() 
  {;}

  void compute(int wid = 0);

 private:
  SinkChannel ch_;
};

class SamsReceive : public kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsReceive(MPILink* link): 
    ch_(link),
    kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH>() 
  {;}

  void compute();

 private:
  SinkChannel ch_;
};

#endif

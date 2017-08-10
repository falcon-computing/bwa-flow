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

// mutex for serializing MPI calls
extern boost::mutex mpi_mutex;

class SeqsDispatch : public kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsDispatch(): kestrelFlow::SinkStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute(int wid = 0);
};

class SeqsReceive : public kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH> {
 public:
  SeqsReceive(): kestrelFlow::SourceStage<SeqsRecord, INPUT_DEPTH>() {;}
  void compute();
};

class SamsSend : public kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsSend(): kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>() {;}
  void compute(int wid = 0);
  std::string serialize(SeqsRecord* data);
};

class SamsReceive : public kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH> {
 public:
  SamsReceive(): kestrelFlow::SourceStage<SeqsRecord, OUTPUT_DEPTH>() {;}
  SeqsRecord deserialize(const char* data, size_t length);
  void compute();
};

#endif

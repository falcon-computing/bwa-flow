#ifndef BWA_FLOW_MPIPIPELINE_H
#define BWA_FLOW_MPIPIPELINE_H

#include "Pipeline.h"
#include "MPIChannel.h"
#include "util.h"

template<typename Record, typename C>
class SendRecord : public kestrelFlow::SinkStage<Record, INPUT_DEPTH> {
 public:
  SendRecord(MPILink* link): 
    ch_(link),
    kestrelFlow::SinkStage<Record, INPUT_DEPTH>() 
  {;}

  void compute(int wid = 0) {
    while (!ch_.sendFinished()) { 
      Record input;
      bool ready = this->getInput(input);

      while (!this->isFinal() && !ready) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
        ready = this->getInput(input);
      }
      if (ready) { // record is a valid new input
        uint64_t start_ts = getUs();

        // Serialize output record
        std::string ser_data = serialize(input);
        int length = ser_data.length();

        freeRecord(input);

        DLOG_IF(INFO, VLOG_IS_ON(2)) << "Serializing one " << input.name 
          << " of " << length / 1024  << "kb"
          << " in " << getUs() - start_ts << " us";

        start_ts = getUs();

        // dispatch data to slaves
        ch_.send(ser_data.c_str(), length);

        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending " << input.name << "-"
          << input.start_idx
          << " takes " << getUs() - start_ts << " us";
      }
      else {
        // this means isFinal() is true and input queue is empty
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finish sending " << input.name
          << ", start sending finish signals";
        ch_.retire();
      }
    }
  }
 private:
  C ch_;
};

template<typename Record, typename C>
class RecvRecord : public kestrelFlow::SourceStage<Record, INPUT_DEPTH> {
 public:
  RecvRecord(MPILink* link): 
    ch_(link),
    kestrelFlow::SourceStage<Record, INPUT_DEPTH>() 
  {;}

  void compute() {
    while (!ch_.recvFinished()) {
      uint64_t start_ts = getUs();

      // request new data from master
      int length = 0;
      char* ser_data = (char*)ch_.recv(length);

      if (length > 0) {
        Record output;
        deserialize(ser_data, length, output);
        free(ser_data);

        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Receive " << output.name 
          << "-" << output.start_idx 
          << " in " << getUs() - start_ts << " us";

        this->pushOutput(output);
      }
    }
  }

 private:
  C ch_;
};


class SeqsDispatch : public SendRecord<SeqsRecord, SourceChannel> {
 public:
  SeqsDispatch(MPILink *link): SendRecord<SeqsRecord, SourceChannel>(link)
  {;}
};

class SeqsGather : public RecvRecord<SeqsRecord, SourceChannel> {
 public:
  SeqsGather(MPILink *link): RecvRecord<SeqsRecord, SourceChannel>(link)
  {;}
};



#endif

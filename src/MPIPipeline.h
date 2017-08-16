#ifndef BWA_FLOW_MPIPIPELINE_H
#define BWA_FLOW_MPIPIPELINE_H

#include "Pipeline.h"
#include "MPIChannel.h"
#include "util.h"

template<typename Record>
class SendStage : public kestrelFlow::SinkStage<Record, INPUT_DEPTH> {
 public:
  SendStage(Channel* ch): 
    ch_(ch), 
    kestrelFlow::SinkStage<Record, INPUT_DEPTH>() 
  {;}

  void compute(int wid = 0) {
    Record r;
    while (!ch_->sendFinished()) { 
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
        ch_->send(ser_data.c_str(), length);

        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending " << input.name << "-"
          << input.start_idx
          << " takes " << getUs() - start_ts << " us";
      }
      else {
        // this means isFinal() is true and input queue is empty
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finish sending " << input.name
          << ", start sending finish signals";
        ch_->retire();
      }
    }
  }
 private:
  Channel* ch_;
};

template<typename Record>
class RecvStage : public kestrelFlow::SourceStage<Record, INPUT_DEPTH> {
 public:
  RecvStage(Channel* ch): 
    ch_(ch),
    kestrelFlow::SourceStage<Record, INPUT_DEPTH>() 
  {;}

  void compute() {
    Record r;
    while (!ch_->recvFinished()) {
      uint64_t start_ts = getUs();

      // request new data from master
      int length = 0;
      char* ser_data = (char*)ch_->recv(length);

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
  Channel* ch_;
};
#endif

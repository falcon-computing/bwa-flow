#include <boost/asio.hpp>
#include <boost/function_types/result_type.hpp>
#include <boost/make_shared.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread.hpp>
#include <fstream>
#include <list>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "mpi.h"

#include "bwa/utils.h"
#include "kflow/Queue.h"

#ifdef USE_HTSLIB
#include "htslib/ksort.h"
#endif

#include "bwa_wrapper.h"
#include "config.h"
#include "MPIPipeline.h"  
#include "util.h"

void SeqsDispatch::compute(int wid) {

  uint64_t num_seqs_produced = 0;

  while (!ch_.isFinished()) { 
    SeqsRecord record;
    bool ready = this->getInput(record);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(record);
    }
    if (ready) { // record is a valid new input
      uint64_t start_ts = getUs();

      // Serialize output record
      std::string ser_data = serialize(record);
      int length = ser_data.length();

      if (length <= 0) {
        throw std::runtime_error("Possible overflow of msg length");
      }

      for (int i = 0; i < record.batch_num; i++) {
        free(record.seqs[i].name);
        free(record.seqs[i].comment);
        free(record.seqs[i].seq);
        free(record.seqs[i].qual);
      }
      free(record.seqs);

      DLOG_IF(INFO, VLOG_IS_ON(2)) << "Serializing seq batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      // dispatch data to slaves
      ch_.send(ser_data.c_str(), length);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending seqs batch " << record.start_idx
        << " takes " << getUs() - start_ts << " us";
    }
    else {
      // this means isFinal() is true and input queue is empty
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Finish reading seqs, start send finish signals";
      ch_.finish();
    }
  }
}

void SeqsReceive::compute() {

  while (!ch_.isFinished()) {
    uint64_t start_ts = getUs();

    // request new data from master
    int length = 0;
    char* ser_data = (char*)ch_.recv(length);

    if (length > 0) {
      SeqsRecord output;
      deserialize(ser_data, length, output);
      free(ser_data);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Receive one read batch in "
        << getUs() - start_ts << " us";

      pushOutput(output);
    }
  }
}

void ChainsToRegionsOffload::compute(int wid) {
  int rank   = 0; //mpi_rank;
  int nprocs = 0; //mpi_nprocs;
  bool finished = false;
  while (!finished) {
    ChainsRecord input;
    bool ready = this->getInput(input);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }

    if (!ready) { 
      // This means isFinal() is true and input queue is empty
      int length = 0;
    }
  }
}

void SamsSend::compute(int wid) {

#ifndef USE_HTSLIB    
  int rank   = 0; // mpi_rank;
  int nprocs = 0; // mpi_nprocs;

  while (!ch_->isFinished()) {
    SeqsRecord input;
    bool ready = this->getInput(input);

    while (!this->isFinal() && !ready) {
      boost::this_thread::sleep_for(boost::chrono::microseconds(100));
      ready = this->getInput(input);
    }

    if (!ready) { 
      // This means isFinal() is true and input queue is empty
      ch_->finish(); 
    }
    else {
      uint64_t start_ts = getUs();

      // Serialize data and send to master
      std::string ser_data = serialize(input);
      int length = ser_data.length();

      if (length <= 0) {
        throw std::runtime_error("Possible overflow of msg length");
      }
      DLOG_IF(INFO, VLOG_IS_ON(2)) << "Serializing sam batch of "
        << length / 1024  << "kb"
        << " in " << getUs() - start_ts << " us";

      start_ts = getUs();

      ch_->send(ser_data.c_str(), length);

      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Sending sam batch " << input.start_idx
        << " to master takes " << getUs() - start_ts << " us";

      //freeSeqs(input.seqs, input.batch_num);
      for (int i = 0; i < input.batch_num; i++) {
        free(input.seqs[i].sam);
      }
      free(input.seqs);
    }
  }
#else
  LOG(ERROR) << "SamSend() is not supported in sorted bwa version";
#endif
}

void SamsReceive::compute() {
#ifndef USE_HTSLIB    
  while (!unfinished_proc.empty()) {

    int length = 0;
    char* ser_data = ch_->recv(length);

    if (length > 0) {
      // Allocate buffer for serialized obj
      SeqsRecord output;
      deserialize(ser_data, length, output);
      free(ser_data);

      pushOutput(output);
    }
  }
#else
  LOG(ERROR) << "SamReceive() is not supported in sorted bwa version";
#endif
}

#include <iostream>
#include <vector>
#include <glog/logging.h>
#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"
#include "kflow/MapPartitionStage.h"

using namespace kestrelFlow;

class RandGenStage : 
  public MapStage<int, std::vector<double>*>
{
public:
  RandGenStage(): MapStage<int, std::vector<double>*>() {;}

  std::vector<double>* compute(int const & l) {

    try {
      boost::any constant = this->getConst("length");
      int length = boost::any_cast<int>(constant);

      std::vector<double>* output = new std::vector<double>(length);

      for (int i=0; i<length; i++) {
        (*output)[i] = (double)i;
      }

      VLOG(1) << "Generated one vector";
      return output;
    }
    catch (paramError &e) {
      LOG(ERROR) << e.what();
      return NULL;
    }
    catch (boost::bad_any_cast &) {
      DLOG(ERROR) << "type of length mismatch";
      return NULL;
    }
  }
};

class NormStage :
  public MapPartitionStage<std::vector<double>*, double> 
{
public:
  NormStage(bool is_dyn, int n):
    MapPartitionStage<std::vector<double>*, double>(n, is_dyn) {;}

  void compute(int wid) {
    if (this->isDynamic()) {
      std::vector<double>* input;
      bool ready = this->getInput(input);

      int counter = 0;
      while (!this->isFinal() && !ready && counter < 20) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
        ready = this->getInput(input);
        counter ++;
      }
      if (!ready) { 
        // this means isFinal() is true or timeout
        return; 
      }

      double norm = 0;
      for (int i=0; i<input->size(); i++) {
        double val = (*input)[i];
        norm += val*val;
      }
      VLOG(1) << "Computed one vector";

      this->pushOutput(norm);
    }
    else {
      while (true) {
        std::vector<double>* input;
        bool ready = this->getInput(input);

        while (!this->isFinal() && !ready) {
          boost::this_thread::sleep_for(boost::chrono::microseconds(100));
          ready = this->getInput(input);
        }
        if (!ready) { 
          // this means isFinal() is true or timeout
          return; 
        }

        double norm = 0;
        for (int i=0; i<input->size(); i++) {
          double val = (*input)[i];
          norm += val*val;
        }
        VLOG(1) << "Computed one vector";

        this->pushOutput(norm);
      }
    }
  }
};


int main(int argc, char** argv) {

  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

  int n = 8;
  int length = 8;
  int stage1_workers = 4;
  int stage2_workers = 1;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    length = atoi(argv[2]);
  }
  if (argc > 3) {
    stage1_workers = atoi(argv[3]); 
  }
  if (argc > 4) {
    stage2_workers = atoi(argv[4]); 
  }

  Pipeline norm_pipeline(2, 2);

  norm_pipeline.addConst("length", length);

  RandGenStage stage1;
  NormStage stage2(false, 1);

  norm_pipeline.addStage(0, &stage1);
  norm_pipeline.addStage(1, &stage2);

  norm_pipeline.start();

  Queue<int>* input_queue = static_cast<Queue<int>*>(
                              norm_pipeline.getInputQueue().get());
  Queue<double>* output_queue = static_cast<Queue<double>*>(
                              norm_pipeline.getOutputQueue().get());

  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 32; k++) {
      input_queue->push(0);
    }
    for (int k = 0; k < 32; k++) {
      double out;
      output_queue->pop(out);

      double correct = 0;
      for (int i = 0; i < length; i++) {
        correct += (double)i*i;
      }
      if (correct != out) {
        LOG(ERROR) << "Results does not match: "
          << correct << " != " << out;
      }
    }
  }
  norm_pipeline.finalize();
  DLOG(INFO) << "Pipeline is finalized";

  // gracefully end the pipeline
  norm_pipeline.wait();

  return 0;
}

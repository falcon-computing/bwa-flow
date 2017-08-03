#include <iostream>
#include <vector>
#include <glog/logging.h>
#include "kflow/Pipeline.h"
#include "kflow/MapStage.h"

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
  public MapStage<std::vector<double>*, double> 
{
public:
  NormStage(): MapStage<std::vector<double>*, double>() {;}

  double compute(std::vector<double>* const & input) {

    double norm = 0;
    for (int i=0; i<input->size(); i++) {
      double val = (*input)[i];
      norm += val*val;
    }
    VLOG(1) << "Computed one vector";
    return norm;
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
  NormStage stage2;

  norm_pipeline.addStage(0, &stage1);
  norm_pipeline.addStage(1, &stage2);

  norm_pipeline.start();

  Queue<int>* input_queue = static_cast<Queue<int>*>(
                              norm_pipeline.getInputQueue().get());
  Queue<double>* output_queue = static_cast<Queue<double>*>(
                              norm_pipeline.getOutputQueue().get());

  for (int i=0; i<n; i+=64) {
    for (int k = i; k < 64 && k < n; k++) {
      input_queue->push(0);
    }
    for (int k = i; k < 64 && k < n; k++) {
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

  // gracefully end the pipeline
  norm_pipeline.wait();

  return 0;
}

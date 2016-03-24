#include <iostream>
#include <vector>
#include "kflow.h"

using namespace kestrelFlow;

class RandGenStage : 
  public MapStage<int, std::vector<double>*>
{
public:
  RandGenStage(int n): MapStage<int, std::vector<double>*>(n) {;}

  std::vector<double>* compute(int const & l) {

    try {
      boost::any constant = this->getConst("length");
      int length = boost::any_cast<int>(constant);

      std::vector<double>* output = new std::vector<double>(length);

      for (int i=0; i<length; i++) {
        (*output)[i] = (double)rand()/RAND_MAX;
      }
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
  NormStage(int n): MapStage<std::vector<double>*, double>(n) {;}

  double compute(std::vector<double>* const & input) {

    double norm = 0;
    for (int i=0; i<input->size(); i++) {
      double val = (*input)[i];
      norm += val*val;
    }
    return norm;
  }
};


int main(int argc, char** argv) {

  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

  int n = 8;
  int length = 8;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    length = atoi(argv[2]);
  }

  Pipeline norm_pipeline(2);

  norm_pipeline.addConst("length", length);

  RandGenStage stage1(8);
  NormStage stage2(8);

  norm_pipeline.addStage(0, &stage1);
  norm_pipeline.addStage(1, &stage2);
  norm_pipeline.start();

  Queue<int>* input_queue = static_cast<Queue<int>*>(
                              norm_pipeline.getInputQueue());
  Queue<double>* output_queue = static_cast<Queue<double>*>(
                              norm_pipeline.getOutputQueue());

  for (int i=0; i<n; i++) {
    input_queue->push(0);
  }
  norm_pipeline.finalize();

  for (int i=0; i<n; i++) {
    double out;
    output_queue->pop(out);
    std::cout << out << std::endl;
  }

  // gracefully end the pipeline
  norm_pipeline.wait();
  
  return 0;
}

#include <iostream>
#include <vector>
#include "kflow.h"

using namespace kestrelFlow;

class RandGenStage : 
  public Stage<int, std::vector<double>*>
{
public:
  RandGenStage(int n): Stage<int, std::vector<double>*>(n) {;}

  std::vector<double>* compute(int const & length) {

    std::vector<double>* output = new std::vector<double>(length);

    for (int i=0; i<length; i++) {
      (*output)[i] = (double)rand()/RAND_MAX;
    }

    return output;
  }
};

class NormStage :
  public Stage<std::vector<double>*, double> 
{
public:
  NormStage(int n): Stage<std::vector<double>*, double>(n) {;}

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

  int n = 8;
  int length = 8;

  if (argc > 1) {
    n = atoi(argv[1]);
  }
  if (argc > 2) {
    length = atoi(argv[2]);
  }

  Pipeline norm_pipeline(2);

  RandGenStage stage1(4);
  NormStage stage2(2);

  norm_pipeline.addStage(0, &stage1);
  norm_pipeline.addStage(1, &stage2);
  norm_pipeline.start();

  Queue<int>* input_queue = static_cast<Queue<int>*>(
                              norm_pipeline.getInputQueue());
  Queue<double>* output_queue = static_cast<Queue<double>*>(
                              norm_pipeline.getOutputQueue());

  for (int i=0; i<n; i++) {
    input_queue->push(length);
  }
  for (int i=0; i<n; i++) {
    double out;
    output_queue->pop(out);
    std::cout << out << std::endl;
  }
  norm_pipeline.stop();
  
  return 0;
}

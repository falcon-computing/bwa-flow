#include <iostream>
#include <vector>
#include "kflow.h"

using namespace kestrelFlow;

template <typename U>
class vector {
public:
  vector(int n): size(n) {
    if (n>0) {
      data = new U[n];
    }
  }
  ~vector() {
    if (data) delete [] data;
  }
  int size;
  U* data;
};

typedef vector<double>* vector_ptr;

class RandGenStage : 
  public Stage<int, 32, vector_ptr, 32>
{
public:
  RandGenStage(int n): Stage<int, 32, vector_ptr, 32>(n) {;}

  void compute() {

    int length = readInput();

    vector_ptr output = new vector<double>(length);

    for (int i=0; i<length; i++) {
      output->data[i] = (double)rand()/RAND_MAX;
    }

    writeOutput(output);
  }
};

class NormStage :
  public Stage<vector_ptr, 32, double, 32> 
{
public:
  NormStage(int n): Stage<vector_ptr, 32, double, 32>(n) {;}

  void compute() {

    vector_ptr input = readInput();
    
    double norm = 0;
    for (int i=0; i<input->size; i++) {
      double val = input->data[i];
      norm += val*val;
    }

    delete input;

    writeOutput(norm);
  }
};


int main() {

  int length = 8;
  int n = 8;

  Pipeline norm_pipeline(2);

  RandGenStage stage1(4);
  NormStage stage2(2);

  norm_pipeline.addStage(0, &stage1);
  norm_pipeline.addStage(1, &stage2);
  norm_pipeline.start();

  Queue<int, 32>* input_queue = dynamic_cast<Queue<int, 32>*>(norm_pipeline.getInputQueue());
  Queue<double, 32>* output_queue = dynamic_cast<Queue<double, 32>*>(norm_pipeline.getOutputQueue());

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

#include <iostream>
#include <fstream>
#include <vector>

#include "kflow/Pipeline.h"
#include "kflow/SourceStage.h"
#include "kflow/SinkStage.h"

#include "mkl.h"
#include "mkl_lapacke.h"

using namespace kestrelFlow;

void mmul(
    double* a,
    double* b,
    double* c,
    int m, int k, int n)
{
  for (int t=0; t<k; t++) {
    for (int j=0; j<n; j++) {
      for (int i=0; i<m; i++) {
        c[i*n+j] += a[i*k+t] * b[t*n+j];
      }
    }
  }
}

struct Record {
  Record(int _m, int _k, int _n):
    m(_m), k(_k), n(_n)
  {
    a = (double*)mkl_malloc(m*k*sizeof(double), 64);
    b = (double*)mkl_malloc(k*n*sizeof(double), 64);
    c = (double*)mkl_malloc(m*n*sizeof(double), 64);
  }
  ~Record() {
    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
  }

  MKL_INT m;
  MKL_INT k;
  MKL_INT n;
  double* a;
  double* b;
  double* c;
};

class Producer : 
  public SourceStage<Record*>
{
public:
  Producer(): SourceStage<Record*>() {;}

  void compute() {

    try {
      boost::any var = this->getConst("m");
      int m = boost::any_cast<int>(var);

      var = this->getConst("k");
      int k = boost::any_cast<int>(var);

      var = this->getConst("n");
      int n = boost::any_cast<int>(var);

      var = this->getConst("batch");
      int batch = boost::any_cast<int>(var);
      
      for (int b=0; b<batch; b++) {
        Record* output = new Record(m, k, n);

        for (int i=0; i<m*k; i++) {
          output->a[i] = (double)i/m*k;
        }
        for (int i=0; i<k*n; i++) {
          output->b[i] = (double)i/n*k;
        }
        for (int i=0; i<m*n; i++) {
          output->c[i] = 0;
        }

        pushOutput(output);
      }
    }
    catch (paramError &e) {
      LOG(ERROR) << e.what();
      return;
    }
    catch (boost::bad_any_cast &) {
      LOG(ERROR) << "type of input_fname mismatch";
      return;
    }
  }
};

class Consumer :
  public SinkStage<Record*>
{
public:
  Consumer(): SinkStage<Record*>() {}

  void compute() {

    while (true) 
    {
      Record* input;
      bool ready = this->getInput(input);

      while (!this->isFinal() && !ready) {
        boost::this_thread::sleep_for(boost::chrono::microseconds(100));
        ready = this->getInput(input);
      }
      if (!ready) { 
        // this means isFinal() is true and input queue is empty
        break; 
      }

      double alpha = 1.0;
      double beta = 0.0;

      MKL_INT m = input->m;
      MKL_INT n = input->n;
      MKL_INT k = input->k;
      double* a = input->a;
      double* b = input->b;
      double* c = input->c;

      uint64_t start_ts = getUs();
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
          m, n, k, alpha, 
			    a, k, 
			    b, n, beta,
			    c, n);
      //mmul(a, b, c, m, k, n);
      LOG(INFO) << "Compute time is " <<  getUs() - start_ts << " us";
    }
  }
};

int main(int argc, char** argv) {

  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

  int batch = 4;

  int m = 1024;
  int k = 1024;
  int n = 1024;

  if (argc > 3) {
    m = atoi(argv[1]);
    k = atoi(argv[2]);
    n = atoi(argv[3]);
  }
  else if (argc > 2) {
    m = atoi(argv[1]);
    k = atoi(argv[1]);
    n = atoi(argv[2]);
  }
  else if (argc > 1) {
    m = atoi(argv[1]);
    k = atoi(argv[1]);
    n = atoi(argv[1]);
  }

  Pipeline mmul(2, 0);

  mmul.addConst("m", m);
  mmul.addConst("k", k);
  mmul.addConst("n", n);
  mmul.addConst("batch", batch);

  Producer p_stage;
  Consumer c_stage;

  mmul.addStage(0, &p_stage);
  mmul.addStage(1, &c_stage);
  mmul.start();

  mmul.wait();
  //mmul.printPerf();

  return 0;
}

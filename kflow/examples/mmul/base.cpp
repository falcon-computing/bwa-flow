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
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      for (int t=0; t<k; t++) {
        c[i*n+j] += a[i*k+t] * b[t*n+j];
      }
    }
  }
}

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

  uint64_t total_t = 0;
  uint64_t gen_t = 0;
  uint64_t compute_t = 0;

  uint64_t start_ts = getUs();
  for (int k=0; k<batch; k++) {

    uint64_t start_ts = getUs();
    double alpha = 1.0;
    double beta = 0.0;

    double* a = (double*)mkl_malloc(m*k*sizeof(double), 64);
    double* b = (double*)mkl_malloc(k*n*sizeof(double), 64);
    double* c = (double*)mkl_malloc(m*n*sizeof(double), 64);

    for (int i=0; i<m*k; i++) {
      a[i] = (double)i/m*k;
    }
    for (int i=0; i<k*n; i++) {
      b[i] = (double)i/n*k;
    }
    for (int i=0; i<m*n; i++) {
      c[i] = 0;
    }

    gen_t += getUs() - start_ts;

    start_ts = getUs();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        m, n, k, alpha, 
        a, k, 
        b, n, beta,
        c, n);
    //mmul(a, b, c, m, k, n);

    LOG(INFO) << "Compute time is " <<  getUs() - start_ts << " us";

    mkl_free(a);
    mkl_free(b);
    mkl_free(c);
  }
  total_t = getUs() - start_ts;

  return 0;
}



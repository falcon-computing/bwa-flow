#include <glog/logging.h>
#include <gtest/gtest.h>

#include "FPGAAgent.h"
#include "TestCommon.h"

gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;
ktp_aux_t aux;

FPGAAgent* agent;

boost::mutex mpi_mutex;

int main(int argc, char *argv[]) {

  google::InitGoogleLogging(argv[0]);
  ::testing::InitGoogleTest(&argc, argv);

  memset(&aux, 0, sizeof(ktp_aux_t));

  // Get the index and the options
  pre_process(argc-1, argv+1, &aux, true);

  int ret = RUN_ALL_TESTS();

  // Free all global variables
  free(aux.opt);
  bwa_idx_destroy(aux.idx);
  kseq_destroy(aux.ks);
  err_gzclose(fp_idx); 
  kclose(ko_read1);

  if (aux.ks2) {
    kseq_destroy(aux.ks2);
    err_gzclose(fp2_read2); kclose(ko_read2);
  }

  return ret;
}



#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <glog/logging.h>
#include <string>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "kflow/Pipeline.h"

#include "bwa_wrapper.h"
#include "Pipeline.h"
#include "util.h"
#include "FPGAAgent.h"

FPGAAgent* agent;

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

int main(int argc, char *argv[]) {

  // Initialize Google Log
  FLAGS_logtostderr = 1;
  google::InitGoogleLogging(argv[0]);

	double t_real = realtime();

  // Preprocessing
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2;
  extern void *ko_read1, *ko_read2;
  ktp_aux_t aux;
  memset(&aux, 0, sizeof(ktp_aux_t));

  kstring_t pg = {0,0,0};
  ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
  for (int i = 1; i < argc; i++) {
    ksprintf(&pg, " %s", argv[i]);
  }
  bwa_pg = pg.s;

  int chunk_size = 2000;

  // Start FPGA agent
  agent = new FPGAAgent("/curr/diwu/prog/acc_lib/bwa-sm/sm-80pe.xclbin", 2000);

  // Get the index and the options
  pre_process(argc-1, argv+1, &aux);

  kestrelFlow::Pipeline bwa_flow(5);

  SeqsProducer    input_stage;
  SeqsToChains    seq2chain_stage(6);
  ChainsToRegions chain2reg_stage(1);
  RegionsToSam    reg2sam_stage(2);
  PrintSam        output_stage;

  bwa_flow.addConst("aux", &aux);
  bwa_flow.addConst("chunk_size", chunk_size);

  bwa_flow.addStage(0, &input_stage);
  bwa_flow.addStage(1, &seq2chain_stage);
  bwa_flow.addStage(2, &chain2reg_stage);
  bwa_flow.addStage(3, &reg2sam_stage);
  bwa_flow.addStage(4, &output_stage);
  bwa_flow.start();
  bwa_flow.wait();
  bwa_flow.printPerf();
  
  // Free all global variables
  delete agent;

  free(aux.opt);
  bwa_idx_destroy(aux.idx);
  kseq_destroy(aux.ks);
  err_gzclose(fp_idx); 
  kclose(ko_read1);

  if (aux.ks2) {
    kseq_destroy(aux.ks2);
    err_gzclose(fp2_read2); kclose(ko_read2);
  }

	err_fflush(stdout);
	err_fclose(stdout);

  fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
  fprintf(stderr, "[%s] CMD:", __func__);
  for (int i = 0; i < argc; ++i) {
    fprintf(stderr, " %s", argv[i]);
  }
  fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());

	free(bwa_pg);
  return 0;
}

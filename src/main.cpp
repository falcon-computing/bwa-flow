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
#include "mpi.h"
#include "kflow/Pipeline.h"

#include "bwa_wrapper.h"
#include "config.h"
#include "FPGAAgent.h"
#include "Pipeline.h"
#include "util.h"

FPGAAgent* agent;

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

int main(int argc, char *argv[]) {

  // Initialize MPI
  int init_ret = MPI::Init_thread(MPI_THREAD_MULTIPLE);

  if (init_ret != MPI_THREAD_MULTIPLE) {
    LOG(ERROR) << "Available thread level is " << init_ret;
    throw std::runtime_error("Cannot initialize MPI with threads");
  }
  int rank   = MPI::COMM_WORLD.Get_rank();
  int nprocs = MPI::COMM_WORLD.Get_size();

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

  int chunk_size = CHUNK_SIZE;

  // Get the index and the options
  pre_process(argc-1, argv+1, &aux, rank==0);

  kestrelFlow::Pipeline bwa_io(2);
  kestrelFlow::Pipeline bwa_flow(5);

  // stages for bwa file in/out
  SeqsDispatcher input_stage;
  PrintSam       output_stage;

  // stages for bwa computation
  SeqsReceiver    recv_stage;
  SeqsToChains    seq2chain_stage(STAGE_1_WORKER_NUM);
  ChainsToRegions chain2reg_stage(STAGE_2_WORKER_NUM);
  RegionsToSam    reg2sam_stage(STAGE_3_WORKER_NUM);
  SendSam         send_stage;

  if (rank == 0) { 
    // Master process start a separate pipeline just to
    // distribute the input seqs and gather output sam
    bwa_io.addConst("aux", &aux);
    bwa_io.addConst("mpi_rank", rank);
    bwa_io.addConst("mpi_nprocs", nprocs);

    bwa_io.addStage(0, &input_stage);
    bwa_io.addStage(1, &output_stage);

    bwa_io.start();
  }
  // Start FPGA agent
  //agent = new FPGAAgent(FPGA_PATH, chunk_size);

  bwa_flow.addConst("aux", &aux);
  bwa_flow.addConst("chunk_size", chunk_size);
  bwa_flow.addConst("mpi_rank", rank);
  bwa_flow.addConst("mpi_nprocs", nprocs);

  bwa_flow.addStage(0, &recv_stage);
  bwa_flow.addStage(1, &seq2chain_stage);
  bwa_flow.addStage(2, &chain2reg_stage);
  bwa_flow.addStage(3, &reg2sam_stage);
  bwa_flow.addStage(4, &send_stage);
  bwa_flow.start();
  bwa_flow.wait();

  // Free all global variables
  //delete agent;

  if (rank == 0) {
    bwa_io.wait();
  }

  free(aux.opt);
  bwa_idx_destroy(aux.idx);

  if (rank == 0) { // master process
    kseq_destroy(aux.ks);
    err_gzclose(fp_idx); 
    kclose(ko_read1);

    if (aux.ks2) {
      kseq_destroy(aux.ks2);
      err_gzclose(fp2_read2); kclose(ko_read2);
    }

    kstring_t pg = {0,0,0};
    ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
    for (int i = 1; i < argc; i++) {
      ksprintf(&pg, " %s", argv[i]);
    }
    bwa_pg = pg.s;

    fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
    fprintf(stderr, "[%s] CMD:", __func__);
    for (int i = 0; i < argc; ++i) {
      fprintf(stderr, " %s", argv[i]);
    }
    fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());

    err_fflush(stdout);
    err_fclose(stdout);

    free(bwa_pg);
  }

  MPI_Finalize();

  return 0;
}

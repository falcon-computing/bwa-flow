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

#ifdef SCALE_OUT
  // Initialize MPI
  int init_ret = MPI::Init_thread(MPI_THREAD_MULTIPLE);

  if (init_ret != MPI_THREAD_MULTIPLE) {
    LOG(ERROR) << "Available thread level is " << init_ret;
    throw std::runtime_error("Cannot initialize MPI with threads");
  }
  int rank   = MPI::COMM_WORLD.Get_rank();
  int nprocs = MPI::COMM_WORLD.Get_size();
#else
  const int rank = 0;
#endif

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

  kestrelFlow::Pipeline scatter_flow(2);
  kestrelFlow::Pipeline gather_flow(2);
  //kestrelFlow::Pipeline compute_flow(5);
  kestrelFlow::Pipeline compute_flow(3);

  // Stages for bwa file in/out
  SeqsRead     input_stage;
  SamsPrint    output_stage;
  SeqsDispatch scatter_stage;
  SamsReceive  gather_stage;

  // Stages for bwa computation
  SeqsReceive     recv_stage;
  SeqsToSams      seq2sam_stage(12);
  //SeqsToChains    seq2chain_stage(STAGE_1_WORKER_NUM);
  //ChainsToRegions chain2reg_stage(STAGE_2_WORKER_NUM);
  //RegionsToSam    reg2sam_stage(STAGE_3_WORKER_NUM);
  SamsSend        send_stage;

  // Bind global vars to each pipeline
  compute_flow.addConst("aux", &aux);
  compute_flow.addConst("chunk_size", chunk_size);

#ifdef SCALE_OUT
  compute_flow.addConst("mpi_rank", rank);
  compute_flow.addConst("mpi_nprocs", nprocs);

  scatter_flow.addConst("aux", &aux);
  scatter_flow.addConst("mpi_rank", rank);
  scatter_flow.addConst("mpi_nprocs", nprocs);

  gather_flow.addConst("aux", &aux);
  gather_flow.addConst("mpi_rank", rank);
  gather_flow.addConst("mpi_nprocs", nprocs);

  scatter_flow.addStage(0, &input_stage);
  scatter_flow.addStage(1, &scatter_stage);

  gather_flow.addStage(0, &gather_stage);
  gather_flow.addStage(1, &output_stage);

  compute_flow.addStage(0, &recv_stage);
  compute_flow.addStage(1, &seq2sam_stage);
  compute_flow.addStage(2, &send_stage);
  //compute_flow.addStage(1, &seq2chain_stage);
  //compute_flow.addStage(2, &chain2reg_stage);
  //compute_flow.addStage(3, &reg2sam_stage);
  //compute_flow.addStage(4, &send_stage);
  
  if (rank == 0) { 
    scatter_flow.start();
    gather_flow.start();
  }
#else
  compute_flow.addStage(0, &input_stage);
  compute_flow.addStage(1, &seq2sam_stage);
  compute_flow.addStage(2, &output_stage);
#endif

  // Start FPGA agent
  //agent = new FPGAAgent(FPGA_PATH, chunk_size);

  compute_flow.start();
  compute_flow.wait();

  // Free all global variables
  //delete agent;

  if (rank == 0) {
#ifdef SCALE_OUT
    scatter_flow.wait();
    gather_flow.wait();
#endif

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

  free(aux.opt);
  bwa_idx_destroy(aux.idx);

#ifdef SCALE_OUT
  MPI_Finalize();
#endif

  return 0;
}

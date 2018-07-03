#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <ctype.h>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <errno.h>
#include <zlib.h>

#include "mpi.h"

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "kflow/Pipeline.h"

#include "bwa_wrapper.h"
#include "config.h"
#include "MPIPipeline.h"
#include "Pipeline.h"
#include "util.h"
#include "allocation_wrapper.h"

#ifndef VERSION
#define VERSION "untracked"
#endif

// use flexlm
#ifdef USELICENSE
#include "falcon-lic/license.h"
#endif

#ifdef BUILD_FPGA
#include "FPGAAgent.h"
#include "FPGAPipeline.h"
BWAOCLEnv* opencl_env;
#endif

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;
ktp_aux_t* aux;

int main(int argc, char *argv[]) {

  // Print arguments for records
  std::stringstream ss;
  for (int i = 0; i < argc; i++) {
    ss << argv[i] << " ";
  }

  std::stringstream version_str;
  version_str << "Falcon BWA-MEM Version: " << VERSION;

  // Initialize Google Flags
  gflags::SetVersionString(version_str.str());
  gflags::SetUsageMessage(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Initialize Google Log
  google::InitGoogleLogging(argv[0]);

  // Initialize MPI
  int init_ret = MPI::Init_thread(MPI_THREAD_SERIALIZED);

  if (init_ret != MPI_THREAD_SERIALIZED) {
    LOG(ERROR) << "Available thread level is " << init_ret;
    throw std::runtime_error("Cannot initialize MPI with threads");
  }

  int rank = MPI::COMM_WORLD.Get_rank();

#ifdef USELICENSE
  namespace fc   = falconlic;
#ifdef DEPLOY_aws
  fc::enable_aws();
#endif
#ifdef DEPLOY_hwc
  fc::enable_hwc();
#endif
  fc::enable_flexlm();

  namespace fclm = falconlic::flexlm;
  fclm::add_feature(fclm::FALCON_DNA);
  int licret = fc::license_verify();
  if (licret != fc::SUCCESS) {
    LOG(ERROR) << "Cannot authorize software usage: " << licret;
    LOG(ERROR) << "Please contact support@falcon-computing.com for details.";
    return -1;
  }
#endif

  // Preprocessing
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2;
  extern void *ko_read1, *ko_read2;
  aux = new ktp_aux_t;
  if (NULL == aux) {
    std::string err_string = "Memory allocation failed";
    if (errno==12)
      err_string += " due to out-of-memory";
    else
      err_string += " due to internal failure"; 
    throw std::runtime_error(err_string);
  }
  memset(aux, 0, sizeof(ktp_aux_t));

  kstring_t pg = {0,0,0};
  ksprintf(&pg, "@PG\tID:falcon-bwa\tPN:bwa\tVN:%s\tCL:%s", VERSION, argv[0]);
  for (int i = 1; i < argc; i++) {
    ksprintf(&pg, " %s", argv[i]);
  }
  bwa_pg = pg.s;

  // Check sanity of input parameters
  int chunk_size = FLAGS_chunk_size;

  // Get output file folder
  std::string sam_dir = FLAGS_output_dir;
  if (!sam_dir.empty()) {
    sam_dir += "/" + boost::asio::ip::host_name()+
      "-" + std::to_string((long long)getpid());

    if (!boost::filesystem::exists(sam_dir)) {
      // Create output folder if it does not exist
      if (!boost::filesystem::create_directories(sam_dir)) {
        LOG(ERROR) << "Cannot create output dir: " << sam_dir;
        return 1;
      }
      if (FLAGS_sort) {
        DLOG_IF(INFO, FLAGS_v >= 1) << "Putting sorted BAM files to " << sam_dir;
      }
      else {
        DLOG_IF(INFO, FLAGS_v >= 1) << "Putting output to " << sam_dir;
      }
    }
  }
  else {
    DLOG_IF(INFO, FLAGS_v >= 1) << "Putting output to stdout";
  }

  // Produce original BWA arguments
  std::vector<const char*> bwa_args;
  if (argc < 2 || strcmp(argv[1], "mem")) {
    LOG(ERROR) << "Expecting 'mem' as the first argument.";
    return 1;
  }
  bwa_args.push_back("mem");

  // Start to pass Google flags through to bwa
  if (FLAGS_v > 1) {
    // Convert glog verbosity to bwa_verbose
    // v< 2 --> bwa_verbose = 2
    // v>=2 --> bwa_verbose = 3
    bwa_args.push_back("-v");
    bwa_args.push_back("3");
  }
  else {
    bwa_args.push_back("-v");
    bwa_args.push_back("1");
  }
  if (FLAGS_M) {
    bwa_args.push_back("-M");
  }
  if (!FLAGS_R.empty()) {
    bwa_args.push_back("-R"); 
    bwa_args.push_back(FLAGS_R.c_str()); 
  }
  if (FLAGS_output_flag < 0 || FLAGS_output_flag > 2) {
    LOG(ERROR) << "Illegal argument of --output_flag";
    return 1;
  }

#ifdef USE_HTSLIB
  bwa_args.push_back("-o");
  bwa_args.push_back(std::to_string((long long)FLAGS_output_flag).c_str()); 
#endif

  // Pass the rest of the record
  for (int i = 2; i < argc; i++) {
    bwa_args.push_back(argv[i]); 
  }

  // If output_dir is set then redirect sam_header to a file
  int stdout_fd;
  if (rank==0 && !sam_dir.empty()) {
    stdout_fd = dup(STDOUT_FILENO);
    std::string fname = sam_dir + "/header";
    freopen(fname.c_str(), "w+", stdout);
  }

  // Parse BWA arguments and generate index and the options
  if (pre_process(bwa_args.size(), (char**)&bwa_args[0], aux, rank==0)) {
    LOG(ERROR) << "Failed to parse BWA arguments";
    return 1;
  }
  free(bwa_pg);

  // broadcast options from master's pre_process
  MPI::COMM_WORLD.Bcast(&aux->opt->flag, 1, MPI::INT, 0);

  // Restore stdout if stdout is redirected
  if (rank==0 && !sam_dir.empty()) {
    fclose(stdout);
    dup2(stdout_fd, STDOUT_FILENO);
    stdout = fdopen(STDOUT_FILENO, "w");
    close(stdout_fd);
  }

  // initialize MPI link
  MPILink link;

  // channels for record send/recv
  // NOTE: the order and number of channels should match 
  // exactly on each processes, because a static id is
  // used to label each communication channels
  SourceChannel seq_ch(&link, 0, true);
  SourceChannel chn_ch(&link, 0, false);
  SinkChannel   reg_ch(&link, 0, false);

  // stages for bwa file in/out
  KseqsRead       kread_stage;
  KseqsToBseqs    k2b_stage;
  SamsReorder     reorder_stage;
  WriteOutput     write_stage(FLAGS_output_nt);

  // stages for bwa computation
  SeqsToChains    seq2chain_stage(FLAGS_stage_1_nt);
  ChainsToRegions chain2reg_stage(FLAGS_stage_2_nt);
  RegionsToSam    reg2sam_stage(FLAGS_stage_3_nt);

  // stages for record communication
  SendStage<SeqsRecord>    seq_send_stage(&seq_ch);
  RecvStage<SeqsRecord>    seq_recv_stage(&seq_ch);
  SendStage<ChainsRecord>  chn_send_stage(&chn_ch);
  RecvStage<ChainsRecord>  chn_recv_stage(&chn_ch);
  SendStage<RegionsRecord> reg_send_stage(&reg_ch);
  RecvStage<RegionsRecord> reg_recv_stage(&reg_ch);

  double t_real = realtime();

  // setup pipelines
#if 0
  kestrelFlow::Pipeline compute_flow(7, FLAGS_t);

  kestrelFlow::Pipeline send_flow(1, 1);
  kestrelFlow::Pipeline recv_flow(1, 1);
  kestrelFlow::Pipeline region_flow(3, FLAGS_t);

  // Bind global vars to each pipeline
  compute_flow.addConst("sam_dir", sam_dir);

  if (rank == 0) { 
    compute_flow.addStage(0, &kread_stage);
    compute_flow.addStage(1, &k2b_stage);
    compute_flow.addStage(2, &seq2chain_stage);
    compute_flow.addStage(3, &chain2reg_stage);
    compute_flow.addStage(4, &reg2sam_stage);
    compute_flow.addStage(5, &reorder_stage);
    compute_flow.addStage(6, &write_stage);

    send_flow.addStage(0, &chn_send_stage);
    recv_flow.addStage(0, &reg_recv_stage);

    compute_flow.diverge(send_flow, 2);
    compute_flow.converge(recv_flow, 4);
  }
  else {
    region_flow.addStage(0, &chn_recv_stage);
    region_flow.addStage(1, &chain2reg_stage);
    region_flow.addStage(2, &reg_send_stage);
  }

  if (rank == 0) {
    compute_flow.start();
    send_flow.start();
    recv_flow.start();
    compute_flow.wait();
  }
  else {
    region_flow.start();
    region_flow.wait();
  }
#else
  kestrelFlow::Pipeline scatter_flow(3, 3);
  kestrelFlow::Pipeline compute_flow(6, FLAGS_t);

  // Bind global vars to each pipeline
  compute_flow.addConst("sam_dir", sam_dir);


  if (rank == 0) {
    scatter_flow.addStage(0, &kread_stage);
    scatter_flow.addStage(1, &k2b_stage);
    scatter_flow.addStage(2, &seq_send_stage);
    
    scatter_flow.start();
  }
  compute_flow.addStage(0, &seq_recv_stage);
  compute_flow.addStage(1, &seq2chain_stage);
  compute_flow.addStage(2, &chain2reg_stage);
  compute_flow.addStage(3, &reg2sam_stage);
  compute_flow.addStage(4, &reorder_stage);
  compute_flow.addStage(5, &write_stage);

  compute_flow.start();
  compute_flow.wait();
#endif
  MPI::COMM_WORLD.Barrier();

  if (rank == 0) {
    std::cerr << "Version: falcon-bwa " << VERSION << std::endl;
    std::cerr << "Real time: " << realtime() - t_real << " sec, "
      << "CPU time: " << cputime() << " sec" 
      << std::endl;

    kseq_destroy(aux->ks);
    err_gzclose(fp_idx); 
    kclose(ko_read1);

    if (aux->ks2) {
      kseq_destroy(aux->ks2);
      err_gzclose(fp2_read2); kclose(ko_read2);
    }
  }
#ifdef USE_HTSLIB
  bam_hdr_destroy(aux->h);
#endif
  free(aux->opt);
  bwa_idx_destroy(aux->idx);
  delete aux;

  MPI_Finalize();

  return 0;
}


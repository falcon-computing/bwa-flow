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
#include <zlib.h>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "kflow/Pipeline.h"

#ifdef SCALE_OUT
#include "mpi.h"
#endif

#include "bwa_wrapper.h"
#include "config.h"
#include "FPGAAgent.h"
#include "Pipeline.h"
#include "util.h"
#include "Extension.h"

OpenCLEnv* opencl_env;

ExtParam** task_batch[5];
bool pend_flag[5];
bool end_flag;
int pend_depth;
boost::mutex mpi_mutex;

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;

int mpi_rank;
int mpi_nprocs;
ktp_aux_t* aux;

// Original BWA parameters
DEFINE_string(R, "", "-R arg in original BWA");

DEFINE_int32(t, boost::thread::hardware_concurrency(),
    "-t arg in original BWA, total number of parallel threads");

DEFINE_bool(M, false, "-M arg in original BWA");

// Parameters
DEFINE_bool(offload, true,
    "Use three compute pipeline stages to enable offloading"
    "workload to accelerators. "
    "If disabled, --use_fpga, --fpga_path will be discard");

DEFINE_bool(use_fpga, false,
    "Enable FPGA accelerator for SmithWaterman computation");

DEFINE_string(fpga_path, "",
    "File path of the SmithWaterman FPGA bitstream");

DEFINE_int32(chunk_size, 2000,
    "Size of each batch send to the FPGA accelerator");

DEFINE_int32(max_fpga_thread, 1,
    "Max number of threads for FPGA worker");

DEFINE_bool(inorder_output, false, 
    "Whether keep the sequential ordering of the sam file");

DEFINE_string(output_dir, "",
    "If not empty the output will be redirect to --output_dir");

DEFINE_int32(stage_1_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 1");

DEFINE_int32(stage_2_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 2");

DEFINE_int32(stage_3_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 3");

int main(int argc, char *argv[]) {

  // Print arguments for records
  std::stringstream ss;
  for (int i = 0; i < argc; i++) {
    ss << argv[i] << " ";
  }

  // Initialize Google Flags
  gflags::SetUsageMessage(argv[0]);
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  // Initialize Google Log
  google::InitGoogleLogging(argv[0]);

#ifdef SCALE_OUT
  // Initialize MPI
  int init_ret = MPI::Init_thread(MPI_THREAD_SERIALIZED);

  if (init_ret != MPI_THREAD_SERIALIZED) {
    LOG(ERROR) << "Available thread level is " << init_ret;
    throw std::runtime_error("Cannot initialize MPI with threads");
  }
  mpi_rank   = MPI::COMM_WORLD.Get_rank();
  mpi_nprocs = MPI::COMM_WORLD.Get_size();

  int rank = mpi_rank;
#else
  const int rank = 0;
#endif

  // Preprocessing
  extern char *bwa_pg;
  extern gzFile fp_idx, fp2_read2;
  extern void *ko_read1, *ko_read2;
  aux = new ktp_aux_t;
  memset(aux, 0, sizeof(ktp_aux_t));

  kstring_t pg = {0,0,0};
  ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
  for (int i = 1; i < argc; i++) {
    ksprintf(&pg, " %s", argv[i]);
  }
  bwa_pg = pg.s;

  // Check sanity of input parameters
  int chunk_size = FLAGS_chunk_size;

  if (FLAGS_offload && FLAGS_use_fpga) {
    VLOG(1) << "Use FPGA in BWA-FLOW";
    boost::filesystem::wpath file_path(FLAGS_fpga_path);
    if (!boost::filesystem::exists(file_path)) {
      LOG(ERROR) << "Cannot find FPGA bitstream at " 
        << FLAGS_fpga_path;
      return 1;
    }
  }
  else {
    FLAGS_use_fpga = false;
  }

  // Get output file folder
  std::string sam_dir = FLAGS_output_dir;
  if (!sam_dir.empty()) {
#ifdef SCALE_OUT
    sam_dir += "/" + boost::asio::ip::host_name()+
      "-" + std::to_string((long long)getpid());
#endif
    // Create output folder if it does not exist
    if (boost::filesystem::create_directories(sam_dir)) {
      VLOG(1) << "Putting sam output to " << sam_dir;
    }
    else {
      LOG(ERROR) << "Cannot create output dir: " << sam_dir;
      return 1;
    }
  }
  else {
    VLOG(1) << "Putting sam output to stdout";
  }

  // Produce original BWA arguments
  std::vector<const char*> bwa_args;
  if (strcmp(argv[1], "mem")) {
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
    bwa_args.push_back("2");
  }
  if (FLAGS_M) {
    bwa_args.push_back("-M");
  }
  if (!FLAGS_R.empty()) {
    bwa_args.push_back("-R"); 
    bwa_args.push_back(FLAGS_R.c_str()); 
  }

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

  // Restore stdout if stdout is redirected
  if (rank==0 && !sam_dir.empty()) {
    fclose(stdout);
    dup2(stdout_fd, STDOUT_FILENO);
    stdout = fdopen(STDOUT_FILENO, "w");
    close(stdout_fd);
  }

	double t_real = realtime();

  int num_compute_stages = 3;
  if (FLAGS_offload) {
    num_compute_stages = 5;
  }

#ifdef SCALE_OUT
  kestrelFlow::Pipeline scatter_flow(2, 0);
  kestrelFlow::Pipeline gather_flow(2, 0);
#endif
  kestrelFlow::Pipeline compute_flow(num_compute_stages, FLAGS_t);

  // Stages for bwa file in/out
  SeqsRead        read_stage;
  SamsPrint       print_stage;
#ifdef SCALE_OUT
  SeqsDispatch    seq_send_stage;
  SeqsReceive     seq_recv_stage;
  SamsSend        sam_send_stage;
  SamsReceive     sam_recv_stage;
#endif
  // Stages for bwa computation
  SeqsToSams      seq2sam_stage(FLAGS_t);
  SeqsToChains    seq2chain_stage(FLAGS_stage_1_nt);
  ChainsToRegions chain2reg_stage(FLAGS_stage_2_nt);
  RegionsToSam    reg2sam_stage(FLAGS_stage_3_nt);

#ifdef SCALE_OUT
  SeqsReceive* input_stage = &seq_recv_stage;
  kestrelFlow::SinkStage<SeqsRecord, OUTPUT_DEPTH>* output_stage;
  if (FLAGS_inorder_output) {
    output_stage = &sam_send_stage;
  }
  else {
    output_stage = &print_stage;
  }
#else
  SeqsRead*  input_stage  = &read_stage;
  SamsPrint* output_stage = &print_stage;
#endif

  // Bind global vars to each pipeline
  compute_flow.addConst("sam_dir", sam_dir);

#ifdef SCALE_OUT
  gather_flow.addConst("sam_dir", sam_dir);

  scatter_flow.addStage(0, &read_stage);
  scatter_flow.addStage(1, &seq_send_stage);

  if (rank == 0) { 
    scatter_flow.start();
    if (FLAGS_inorder_output) {
      gather_flow.addStage(0, &sam_recv_stage);
      gather_flow.addStage(1, &print_stage);
      gather_flow.start();
    }
  }
#endif
  
  compute_flow.addStage(0, input_stage);

  if (FLAGS_offload) {
    compute_flow.addStage(1, &seq2chain_stage);
    compute_flow.addStage(2, &chain2reg_stage);
    compute_flow.addStage(3, &reg2sam_stage);
    compute_flow.addStage(4, output_stage);
  }
  else {
    compute_flow.addStage(1, &seq2sam_stage);
    compute_flow.addStage(2, output_stage);
  }
  
  // Start FPGA context
  
  boost::thread_group driver_threads; 
  FPGAAgent* agent;
  if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
    try {
      opencl_env = new OpenCLEnv(FLAGS_fpga_path.c_str(), "sw_top");
      //agent = new FPGAAgent(FLAGS_fpga_path.c_str(), chunk_size);
      VLOG(1) << "Configured FPGA bitstream from " << FLAGS_fpga_path;
    }
    catch (std::runtime_error &e) {
      LOG(ERROR) << "Cannot configured FPGA bitstream from " << FLAGS_fpga_path
        << " because: " << e.what();
      return 1;
    }
    for (int i = 0; i < 5; i++) {
      task_batch[i] = new ExtParam*[chunk_size];
      pend_flag[i] = false;
    }
    pend_depth = 0;
    end_flag =false;
    // Create FPGAAgent 
    agent = new FPGAAgent(opencl_env, chunk_size);
    driver_threads.create_thread(boost::bind(&fpga_driver,agent)) ;
  }

  
  compute_flow.start();
  compute_flow.wait();

  if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
    driver_threads.interrupt_all();
    driver_threads.join_all();
    delete opencl_env;
    delete agent;
  }

#ifdef SCALE_OUT
  MPI::COMM_WORLD.Barrier();
  if (rank == 0) {
    scatter_flow.wait();
    if (FLAGS_inorder_output) {
      gather_flow.wait();
    }
#else
  if (rank == 0) {
#endif

    kseq_destroy(aux->ks);
    err_gzclose(fp_idx); 
    kclose(ko_read1);

    if (aux->ks2) {
      kseq_destroy(aux->ks2);
      err_gzclose(fp2_read2); kclose(ko_read2);
    }

    kstring_t pg = {0,0,0};
    ksprintf(&pg, "@PG\tID:bwa\tPN:bwa\tVN:%s\tCL:%s", PACKAGE_VERSION, argv[0]);
    for (int i = 1; i < argc; i++) {
      ksprintf(&pg, " %s", argv[i]);
    }
    bwa_pg = pg.s;

    LOG(INFO) << "Version: " << PACKAGE_VERSION;
    LOG(INFO) << "Command: " << ss.str();
    LOG(INFO) << "Real time: " << realtime() - t_real << " sec, "
              << "CPU time: " << cputime() << " sec";

    err_fflush(stdout);
    err_fclose(stdout);

    free(bwa_pg);
  }

  free(aux->opt);
  bwa_idx_destroy(aux->idx);
  delete aux;

#ifdef SCALE_OUT
  MPI_Finalize();
#endif

  return 0;
}

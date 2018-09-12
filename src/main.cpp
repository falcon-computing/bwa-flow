#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <ctype.h>
#include <gflags/gflags.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <errno.h>
#include <zlib.h>

#ifdef NDEBUG
#define LOG_HEADER "falcon-bwa"
#endif
#include <glog/logging.h>

#include "bwa/bntseq.h"
#include "bwa/bwa.h"
//#include "bwa/bwamem.h"
#include "bwa/kseq.h"
#include "bwa/kvec.h"
#include "bwa/utils.h"
#include "kflow/Pipeline.h"
#include "kflow/MegaPipe.h"

#ifndef VERSION
#define VERSION "untracked"
#endif

#include "bwa_wrapper.h"
#include "config.h"
#include "Pipeline.h"
#include "util.h"
#include "allocation_wrapper.h"

#ifdef BUILD_FPGA
#include "FPGAAgent.h"
#include "FPGAPipeline.h"
BWAOCLEnv* opencl_env;
#endif

// use flexlm
#ifdef USELICENSE
#include "falcon-lic/license.h"
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
  google::SetVersionString(version_str.str());
  google::SetUsageMessage(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);

  // Initialize Google Log
  google::InitGoogleLogging(argv[0]);
  DLOG(INFO) << ss.str();

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


  // Pack all the bwa mem args
  pack_bwa_mem_args( bwa_args );
  // Pass the rest of the record
  for (int i = 2; i < argc; i++) {
    bwa_args.push_back(argv[i]); 
  }

  // If output_dir is set then redirect sam_header to a file
  int stdout_fd;
  if (!sam_dir.empty()) {
    stdout_fd = dup(STDOUT_FILENO);
    std::string fname = sam_dir + "/header";
    freopen(fname.c_str(), "w+", stdout);
  }

  // Parse BWA arguments and generate index and the options
  if (pre_process(bwa_args.size(), (char**)&bwa_args[0], aux, true)) {
    LOG(ERROR) << "Failed to parse BWA arguments";
    return 1;
  }

  // Set up OpenCL environment
#ifdef BUILD_FPGA
  if (FLAGS_offload && FLAGS_use_fpga) {
    boost::filesystem::wpath sw_file_path(FLAGS_sw_fpga_path);
    boost::filesystem::wpath smem_file_path(FLAGS_smem_fpga_path);
    
    if (!boost::filesystem::exists(sw_file_path)   && 
        !boost::filesystem::exists(smem_file_path)    ) {
      DLOG(ERROR) << "Cannot find SW FPGA bitstream at " << FLAGS_sw_fpga_path;
      DLOG(ERROR) << "Cannot find SMem FPGA bitstream at " << FLAGS_smem_fpga_path;
      DLOG(WARNING) << "Continue without using FPGA";
      FLAGS_use_fpga = false;
    }
    if (!boost::filesystem::exists(sw_file_path)) FLAGS_sw_fpga_path = "";
    if (!boost::filesystem::exists(smem_file_path)) FLAGS_smem_fpga_path = "";
    try {
      opencl_env = new BWAOCLEnv();
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "Configured FPGA bitstream from - ";
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "  sw:   " << FLAGS_sw_fpga_path; 
      DLOG_IF(INFO, VLOG_IS_ON(1)) << "  smem: " << FLAGS_smem_fpga_path; 
    }
    catch (std::runtime_error &e) {
      DLOG(ERROR) << "Cannot initialize BWA OpenCL environment";
      DLOG(ERROR) << "because: " << e.what();
      DLOG(WARNING) << "Continue without using FPGA";
      FLAGS_use_fpga = false;
    }
  }
  else {
    FLAGS_use_fpga = false;
  }
  if (FLAGS_use_fpga)
    DLOG_IF(INFO, VLOG_IS_ON(1)) << "Use FPGA in BWA-FLOW";
  else
    FLAGS_max_fpga_thread = 0;
#endif

  // dump aux env

  // Restore stdout if stdout is redirected
  if (!sam_dir.empty()) {
    fclose(stdout);
    dup2(stdout_fd, STDOUT_FILENO);
    stdout = fdopen(STDOUT_FILENO, "w");
    close(stdout_fd);
  }

  double t_real = realtime();

  int num_compute_stages = 8;

  int num_threads = FLAGS_t - FLAGS_extra_thread;
#ifdef BUILD_FPGA
  int smem_fpga_thread = (opencl_env)?opencl_env->smem_fpga_thread_:0;
  int sw_fpga_thread = (opencl_env)?opencl_env->sw_fpga_thread_:0;
  if (FLAGS_use_fpga) num_threads -= (smem_fpga_thread + sw_fpga_thread );
#endif
  kestrelFlow::Pipeline compute_flow(num_compute_stages, num_threads);

  DLOG(INFO) << "Using " << num_threads << " threads for cpu";
#ifdef BUILD_FPGA
  DLOG(INFO) << "Using " << smem_fpga_thread + sw_fpga_thread << " for fpga";
#endif

  // Stages for bwa file in/out
  KseqsRead       kread_stage;
  KseqsToBseqs    k2b_stage(FLAGS_t);
  SamsReorder     reorder_stage;
  SamsSort        sort_stage(FLAGS_t);
  WriteOutput     write_stage(FLAGS_output_nt);

  // Stages for bwa computation
  SeqsToChains     seq2chain_stage(FLAGS_stage_1_nt);
  ChainsToRegions  chain2reg_stage(FLAGS_stage_2_nt);
  RegionsToSam     reg2sam_stage(FLAGS_stage_3_nt);
#ifdef BUILD_FPGA
  // Stages for FPGA acceleration of stage_1
  SeqsToChainsFPGA      seq2chain_fpga_stage(smem_fpga_thread);
  // Stages for FPGA acceleration of stage_2
  ChainsToRegionsFPGA   chain2reg_fpga_stage(sw_fpga_thread);
  //SeqsToChainsFPGA      seq2chain_fpga_stage(std::max(opencl_env->smem_fpga_thread_, 1));
  //ChainsToRegionsFPGA   chain2reg_fpga_stage(std::max(opencl_env->sw_fpga_thread_,1));
  seq2chain_fpga_stage.cpu_stage=&seq2chain_stage;
  chain2reg_fpga_stage.cpu_stage=&chain2reg_stage;
#endif

  kestrelFlow::MegaPipe  bwa_flow_pipe(num_threads, FLAGS_max_fpga_thread);

  try {
    // Bind global vars to each pipeline
    compute_flow.addConst("sam_dir", sam_dir);
  
    compute_flow.addStage(0, &kread_stage);
    compute_flow.addStage(1, &k2b_stage);
    compute_flow.addStage(2, &seq2chain_stage);
    //compute_flow.addStage(2, &seq2chain_fpga_stage);
    compute_flow.addStage(3, &chain2reg_stage);
    compute_flow.addStage(4, &reg2sam_stage);
    compute_flow.addStage(5, &reorder_stage);
    compute_flow.addStage(6, &sort_stage);
    compute_flow.addStage(7, &write_stage);

    bwa_flow_pipe.addPipeline(&compute_flow, 1);
  
#ifdef BUILD_FPGA
    if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
      if (smem_fpga_thread > 0)
        compute_flow.addAccxBckStage(2, &seq2chain_fpga_stage, 2);
      if (sw_fpga_thread > 0)
        compute_flow.addAccxBckStage(3, &chain2reg_fpga_stage, 8);
    }
#endif
    
    t_real = realtime();
    bwa_flow_pipe.start();
    bwa_flow_pipe.wait();
  
#ifdef BUILD_FPGA
    if (FLAGS_use_fpga) {
      // Stop FPGA context
      try {
        delete opencl_env;
      }
      catch (std::runtime_error &e) {
        LOG_IF(ERROR, VLOG_IS_ON(1)) << "Failed to run bwa-flow on FPGA";
        DLOG(ERROR) << "because " << e.what();
        throw e;
      }
    }
#endif
  
    kseq_destroy(aux->ks);
    err_gzclose(fp_idx); 
    kclose(ko_read1);
  
    if (aux->ks2) {
      kseq_destroy(aux->ks2);
      err_gzclose(fp2_read2); kclose(ko_read2);
    }
  
    std::cerr << "Version: falcon-bwa " << VERSION << std::endl;
    std::cerr << "Real time: " << realtime() - t_real << " sec, "
              << "CPU time: " << cputime() << " sec" 
              << std::endl;
  
    //err_fflush(stdout);
    //err_fclose(stdout);
  
    free(bwa_pg);
  
#ifdef USE_HTSLIB
    bam_hdr_destroy(aux->h);
#endif
    free(aux->opt);
    bwa_idx_destroy(aux->idx);
    delete aux;
  }
  catch (...) {
    LOG(ERROR) << "Encountered an internal issue.";
    LOG(ERROR) << "Please contact support@falcon-computing.com for details.";
    return 1;
  }

  return 0;
}

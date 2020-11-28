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
#include <execinfo.h>
#include <signal.h>
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

#include "allocation_wrapper.h"
#include "bwa_wrapper.h"
#include "config.h"
#include "Pipeline.h"
#include "util.h"
#include "BucketSortStage.h"
#include "MarkDupStage.h"
#include "MarkDupPartStage.h"
#include "IndexGenStage.h"
#include "BamReadStage.h"
#include "BamSortStage.h"
#include "BamWriteStage.h"
#include "ReorderAndWriteStage.h"

#ifdef BUILD_FPGA
#include "FPGAAgent.h"
#include "FPGAPipeline.h"
BWAOCLEnv* opencl_env;
#endif

// global parameters
gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;
ktp_aux_t* aux;

#ifndef NDEBUG
void trace_dumper(int sig) {
  void *bt_array[64];
  size_t size;
  size = backtrace(bt_array, 64);
  fprintf(stderr, "Debug: caught error signal%d:\n", sig);
  backtrace_symbols_fd(bt_array, size, STDERR_FILENO);
  exit(1);
}
#endif

int main(int argc, char *argv[]) {
#ifndef NDEBUG
  signal(SIGSEGV, trace_dumper);
  fprintf(stderr, "Debug: registered trace dumper for SIGSEGV\n");
#endif

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
  std::string sam_dir = FLAGS_temp_dir;
  if (!sam_dir.empty()) {
    if (!boost::filesystem::exists(sam_dir)) {
      // Create output folder if it does not exist
      if (!boost::filesystem::create_directories(sam_dir)) {
        LOG(ERROR) << "Cannot create temp dir: " << sam_dir;
        return 1;
      }
      DLOG_IF(INFO, FLAGS_v >= 1) << "Putting temp output to " << sam_dir;
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
    boost::filesystem::wpath fpga_file_path(FLAGS_fpga_path);
    
    if (!boost::filesystem::exists(fpga_file_path)) {
      DLOG(ERROR) << "Cannot find FPGA bitstream at " << FLAGS_fpga_path;
      DLOG(WARNING) << "Continue without using FPGA";
      FLAGS_use_fpga = false;
    }
    else {
      try {
        opencl_env = new BWAOCLEnv();
        DLOG_IF(INFO, VLOG_IS_ON(1)) << "Configured FPGA bitstream from "
                                     << FLAGS_fpga_path; 
      }
      catch (std::runtime_error &e) {
        DLOG(ERROR) << "Cannot initialize BWA OpenCL environment";
        DLOG(ERROR) << "because: " << e.what();
        DLOG(WARNING) << "Continue without using FPGA";
        FLAGS_use_fpga = false;
      }
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

  int if_markdup = 0;
  if (!FLAGS_disable_markdup) {
    if_markdup = 1;
  }
  int if_bucketsort = 0;
  if (!FLAGS_disable_bucketsort) {
    if_bucketsort = 1;
  }
  int num_compute_stages = 9 + if_markdup - if_bucketsort;

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
  SeqsToChains      seq2chain_stage(FLAGS_stage_1_nt);
  ChainsPipe        chainpipe_stage(FLAGS_stage_1_nt);
  ChainsToRegions   chain2reg_stage(FLAGS_stage_2_nt);
  RegionsToSam      reg2sam_stage(FLAGS_stage_3_nt);

  // Stages for markduplicates
  //Markdup           md_stage(FLAGS_stage_3_nt, aux);
  MarkDupStage      md_stage(FLAGS_stage_3_nt, aux);
  MarkDupPartStage  md_part_stage(aux);
  BucketSortStage bucketsort_stage(aux, sam_dir, FLAGS_num_buckets, FLAGS_stage_3_nt);

#ifdef BUILD_FPGA
  // Stages for FPGA acceleration of stage_1
  SeqsToChainsFPGA      seq2chain_fpga_stage(smem_fpga_thread, &seq2chain_stage);
  // Stages for FPGA acceleration of stage_2
  ChainsToRegionsFPGA   chain2reg_fpga_stage(sw_fpga_thread, &chain2reg_stage);
#endif
  
  kestrelFlow::MegaPipe  bwa_flow_pipe(num_threads, FLAGS_max_fpga_thread);

  try {
    // Bind global vars to each pipeline
    compute_flow.addConst("sam_dir", sam_dir);
  
    compute_flow.addStage(0, &kread_stage);
    compute_flow.addStage(1, &k2b_stage);
    if (!FLAGS_disable_smem_cpu)
      compute_flow.addStage(2, &seq2chain_stage);
    else
#ifdef BUILD_FPGA
      if (smem_fpga_thread > 0)
        compute_flow.addStage(2, &seq2chain_fpga_stage);
      else {
        LOG(ERROR) << "broken pipeline at seq2chain stage";
        throw std::runtime_error("broken pipeline at seq2chain stage");
      } 
#else
    {
      LOG(ERROR) << "broken pipeline at seq2chain stage";
      throw std::runtime_error("broken pipeline at seq2chain stage");
    }
#endif
    compute_flow.addStage(3, &chainpipe_stage);
    if (!FLAGS_disable_sw_cpu)
      compute_flow.addStage(4, &chain2reg_stage);
    else
#ifdef BUILD_FPGA
      if (sw_fpga_thread > 0)
        compute_flow.addStage(4, &chain2reg_fpga_stage);
      else {
        LOG(ERROR) << "broken pipeline at chain2reg stage";
        throw std::runtime_error("broken pipeline at chain2reg stage");
      }
#else
    {
      LOG(ERROR) << "broken pipeline at chain2reg stage";
      throw std::runtime_error("broken pipeline at chain2reg stage");
    }
#endif
    compute_flow.addStage(5, &reg2sam_stage);

    if (!FLAGS_disable_markdup) {
      if (FLAGS_inorder_output) {
        compute_flow.addStage(6, &reorder_stage);
        compute_flow.addStage(7, &md_part_stage);
      }
      else {
        compute_flow.addStage(6, &md_stage);
        compute_flow.addStage(7, &reorder_stage);
      }
    }
    else {
      compute_flow.addStage(6, &reorder_stage); 
    }
    if (!FLAGS_disable_bucketsort) {
      compute_flow.addStage(7 + if_markdup, &bucketsort_stage);
    }
    else {
      compute_flow.addStage(7 + if_markdup, &sort_stage);
      compute_flow.addStage(8 + if_markdup, &write_stage);
    }
    
    bwa_flow_pipe.addPipeline(&compute_flow, 1);
#ifdef BUILD_FPGA
    if (FLAGS_use_fpga && FLAGS_max_fpga_thread) {
      if (smem_fpga_thread > 0 && !FLAGS_disable_smem_cpu)
        compute_flow.addAccxBckStage(2, &seq2chain_fpga_stage, 2.5);
      if (sw_fpga_thread > 0 && !FLAGS_disable_sw_cpu)
        compute_flow.addAccxBckStage(4, &chain2reg_fpga_stage, 10);
    }
#endif
    
    t_real = realtime();
    bwa_flow_pipe.start();
    bwa_flow_pipe.wait();
  
    bucketsort_stage.closeBuckets();

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

    free(bwa_pg);
    free(aux->opt);
    bwa_idx_destroy(aux->idx);

    std::cerr << "bwa stage time: " 
              << realtime() - t_real 
              << " s" << std::endl;

    t_real = realtime();
    
    if (!FLAGS_disable_bucketsort && FLAGS_merge_bams) {

      kestrelFlow::Pipeline sort_pipeline(4, num_threads);

      // start sorting buckets and merge them
      IndexGenStage     indexgen_stage(
          FLAGS_num_buckets + (!FLAGS_filter_unmap));
      BamReadStage      bamread_stage(sam_dir, aux->h, FLAGS_t);
      BamSortStage      bamsort_stage(FLAGS_t);
      //ReorderAndWriteStage  reorderwrite_stage(FLAGS_output, aux->h);
      BamWriteStage     bamwrite_stage(
          FLAGS_num_buckets + (!FLAGS_filter_unmap),
          sam_dir, FLAGS_output,
          aux->h, FLAGS_t);

      sort_pipeline.addStage(0, &indexgen_stage);
      sort_pipeline.addStage(1, &bamread_stage);
      sort_pipeline.addStage(2, &bamsort_stage);
      //sort_pipeline.addStage(2 + if_sort, &reorderwrite_stage);
      sort_pipeline.addStage(3, &bamwrite_stage);

      kestrelFlow::MegaPipe mp(num_threads, 0);
      mp.addPipeline(&sort_pipeline, 1);

      mp.start();
      mp.wait();

      std::cerr << "sort stage time: " 
        << realtime() - t_real 
        << " s" << std::endl;
    }
#ifdef USE_HTSLIB
    bam_hdr_destroy(aux->h);
#endif
    delete aux;

    std::cerr << "Version: falcon-bwa " << VERSION << std::endl;
  }
  catch (...) {
    LOG(ERROR) << "Encountered an internal issue.";
    LOG(ERROR) << "Please contact support@falcon-computing.com for details.";
    return 1;
  }

  // delete temp_dir
  if (!FLAGS_disable_bucketsort && FLAGS_merge_bams) {
    boost::filesystem::remove_all(sam_dir);
  }

  return 0;
}

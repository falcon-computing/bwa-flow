#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

#include "TestCommon.h"

// Parameters
DEFINE_bool(offload, true,
    "Use three compute pipeline stages to enable offloading"
    "workload to accelerators. "
    "If disabled, --use_fpga, --fpga_path will be discard");

DEFINE_int32(filter, 0, "Filtering out records with INT bit set"
    "on the FLAG field, similar to the -F argument in samtools");

DEFINE_bool(use_fpga, false,
    "Enable FPGA accelerator for SmithWaterman computation");

DEFINE_bool(sort, true,
    "Enable in-memory sorting of output bam file");

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

DEFINE_int32(output_nt, 3,
    "Total number of parallel threads to use for output stage");

DEFINE_int32(output_flag, 1, 
    "Flag to specify output format: "
    "0: BAM (compressed); 1: BAM (uncompressed); 2: SAM");

DEFINE_int32(max_batch_records, 20, 
    "Flag to specify how many batch to buffer before sort");

gzFile fp_idx, fp2_read2 = 0;
void *ko_read1 = 0, *ko_read2 = 0;
ktp_aux_t* aux;

// bseq_t files
bseq1_t*  g_seqs;
int       g_batch_num;

boost::mutex mpi_mutex;

int    g_argc = 0;
char** g_argv = 0;

int main(int argc, char *argv[]) {

#ifdef USE_MPI
  // Initialize MPI
  int init_ret = MPI::Init_thread(MPI_THREAD_SERIALIZED);

  if (init_ret != MPI_THREAD_SERIALIZED) {
    LOG(ERROR) << "Available thread level is " << init_ret;
    throw std::runtime_error("Cannot initialize MPI with threads");
  }
#endif

  google::InitGoogleLogging(argv[0]);
  ::testing::InitGoogleTest(&argc, argv);

  // suppress messages
  bwa_verbose = 0;

  aux = (ktp_aux_t*)malloc(sizeof(ktp_aux_t));
  memset(aux, 0, sizeof(ktp_aux_t));

  // initialize idx
  if (pre_process(argc-1, argv+1, aux, true)) {
    LOG(ERROR) << "Cannot initialize";
  }

  // save seqs for all tests
  g_seqs = bseq_read(aux->actual_chunk_size, &g_batch_num, aux->ks, aux->ks2);
  if (!g_seqs) {
    throw std::runtime_error("cannot read sequence");
  }

  // run all tests
  int ret = RUN_ALL_TESTS();

  // free bseq1_t
  for (int i = 0; i < g_batch_num; i++) {
    free(g_seqs[i].name); 
    free(g_seqs[i].comment);
    free(g_seqs[i].seq); 
    free(g_seqs[i].qual); 
  }
  free(g_seqs);

  // free all global variables
  free(aux->opt);
  bwa_idx_destroy(aux->idx);
  kseq_destroy(aux->ks);
  err_gzclose(fp_idx); 
  kclose(ko_read1);

  if (aux->ks2) {
    kseq_destroy(aux->ks2);
    err_gzclose(fp2_read2); kclose(ko_read2);
  }
#ifdef USE_MPI
  MPI_Finalize();
#endif

  return ret;
}



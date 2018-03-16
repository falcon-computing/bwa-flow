#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>

#include "config.h"

// Original BWA parameters
DEFINE_string(R, "", "-R arg in original BWA");

DEFINE_int32(t, boost::thread::hardware_concurrency(),
    "-t arg in original BWA, total number of parallel threads");

DEFINE_bool(M, true, "-M arg in original BWA");

// Parameters
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

DEFINE_int32(extra_thread, 1,
    "Adjustment to the total threads");

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

DEFINE_int32(output_nt, 2,
    "Total number of parallel threads to use for output stage");

DEFINE_int32(output_flag, 1, 
    "Flag to specify output format: "
    "0: BAM (compressed); 1: BAM (uncompressed); 2: SAM");

DEFINE_int32(max_batch_records, 40, 
    "Flag to specify how many batch to buffer before sort");

// deprecated
DEFINE_bool(offload, true,
    "(deprecated) Use three compute pipeline stages to enable offloading"
    "workload to accelerators. "
    "If disabled, --use_fpga, --fpga_path will be discard");

DEFINE_string(pac_path, "",
    "(deprecated) File path of the modified reference pac file");



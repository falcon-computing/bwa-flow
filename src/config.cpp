#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/thread.hpp>
#include <gflags/gflags.h>
#include <glog/logging.h>

#include "config.h"

// Original BWA parameters
DEFINE_string(R, "@RG\\tID:sample\\tSM:sample\\tPL:illumina\\tLB:sample",
    "-R arg in original BWA");

DEFINE_int32(t, boost::thread::hardware_concurrency(),
    "-t arg in original BWA, total number of parallel threads");

DEFINE_bool(M, true, "-M arg in original BWA");

DEFINE_int32(k, 0, "-k arg in original BWA, minimum seed length");

DEFINE_bool(1, false, "-1 flag in original BWA, flag for *no_mt_io*");

DEFINE_string(x, "",
              "-x arg in original BWA, read type. Setting -x changes multiple parameters unless overriden");

DEFINE_int32(w, 0, 
             "-w arg in original BWA, band width for banded alignment");

DEFINE_int32(A, 0,
             "-A arg in original BWA, score for a sequence match, which scales options -TdBOELU unless overridden");

DEFINE_int32(B, 0,
             "-B arg in original BWA, penalty for a mismatch");

DEFINE_int32(T, 0,
             "-T arg in original BWA, minimum score to output");

DEFINE_int32(U, 0,
             "-U arg in original BWA, penalty for an unpaired read pair");

DEFINE_bool(P, false,
            "-P flag in original BWA, skip pairing; mate rescue performed unless -S also in use");

DEFINE_bool(a, false,
            "-a flag in original BWA, output all alignments for SE or unpaired PE");

DEFINE_bool(p, false,
            "-p flag in original BWA, smart pairing (ignoring in2.fq)");

DEFINE_bool(S, false,
            "-S flag in original BWA, skip mate rescue");

DEFINE_bool(Y, false,
            "-Y flag in original BWA, use soft clipping for supplementary alignments");

DEFINE_bool(V, false,
            "-V flag in original BWA, output the reference FASTA header in the XR tag");

DEFINE_int32(c, 0,
           "-c INT in original BWA, skip seeds with more than INT occurrences");

DEFINE_int32(d, 0,
             "-d arg in original BWA, off-diagonal X-dropoff");

#if 0
//Overrided by FLAGS_v from glog
DEFINE_bool(v, 1, "-v arg in original BWA, verbose level")
#endif

DEFINE_bool(j, false,
            "-j flag in original BWA, treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)");

DEFINE_double(r, 0,
             "-r FLOAT in original BWA, look for internal seeds inside a seed longer than {-k} * FLOAT");

DEFINE_double(D, 0,
             "-D FLOAT in original BWA, drop chains shorter than FLOAT fraction of the longest overlapping chain");

DEFINE_int32(m, 0,
             "-m INT in original BWA, perform at most INT rounds of mate rescues for each read");

DEFINE_int32(s, 0,
             "(deprecated) -s INT in original BWA, look for internal seeds inside a seed with less than INT occ");

DEFINE_int32(G, 0,
             "-G arg in original BWA"); // no description in original bwa codes

DEFINE_int32(N, 0, 
             "-N arg in original BWA"); // no description in original bwa codes

DEFINE_int32(W, 0,
             "-M INT in original BWA, discard a chain if seeded bases shorter than INT");

DEFINE_int32(y, 0,
             "-y arg in original BWA, seed occurrence for the 3rd round seeding");

DEFINE_bool(C, false,
             "-C flag in original BWA, append FASTA/FASTQ comment to SAM output");

DEFINE_int32(K, 0,
             "-K arg in original BWA, fixed chunk size"); // no description in original bwa codes

DEFINE_double(X, 0,
              "-X arg in original BWA"); // no description in original bwa codes

#ifdef USE_HTSLIB
DEFINE_int32(o, 0,
             "-o arg in original BWA, 0 - BAM (compressed), 1 - BAM (uncompressed), 2 - SAM");
#endif

DEFINE_string(h, "",
              "-h arg in original BWA, if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]");

DEFINE_int32(Q, 0,
             "-Q arg in original BWA");

DEFINE_string(O, "",
             "-O arg in original BWA, gap open penalties for deletions and insertions");

DEFINE_string(E, "",
              "-E arg in original BWA, gap extension penalty; a gap of size k cost '{-O} + {-E}*k'");

DEFINE_string(L, "",
              "-L arg in original BWA, penalty for 5'- and 3'-end clipping");

DEFINE_string(H, "",
              "-H STR/FILE in original BWA, insert STR to header if it starts with @; or insert lines in FILE");

DEFINE_string(I, "",
              "-I arg in original BWA\nspecify the mean, standard deviation (10%% of the mean if absent), max\n(4 sigma from the mean if absent) and min of the insert size distribution.\nFR orientation only.");

// Parameters
DEFINE_int32(filter, 0, "Filtering out records with INT bit set"
    "on the FLAG field, similar to the -F argument in samtools");

DEFINE_bool(use_fpga, true,
    "Enable FPGA accelerator for SMem & SmithWaterman computation");

DEFINE_bool(disable_smem_fpga, false,
    "Disable FPGA accelerator for SMem computation");

DEFINE_bool(disable_smem_cpu, false,
    "Disable CPU for SMem computation");

DEFINE_bool(disable_sw_fpga, false,
    "Disable FPGA accelerator for SmithWaterman computation");

DEFINE_bool(disable_sw_cpu, false,
    "Disable CPU for SmithWaterman computation");

DEFINE_bool(sort, true,
    "(deprecated) Enable in-memory sorting of output bam file");

DEFINE_string(fpga_path, "",
    "File path of the SmithWaterman FPGA bitstream");

DEFINE_int32(chunk_size, 2000,
    "Size of each batch send to the FPGA accelerator");

DEFINE_int32(max_fpga_thread, -1,
    "Max number of threads for FPGA worker");

DEFINE_int32(extra_thread, 1,
    "Adjustment to the total threads");

DEFINE_bool(inorder_output, false, 
    "Whether keep the sequential ordering of the sam file");

DEFINE_int32(stage_1_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 1");

DEFINE_int32(stage_2_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 2");

DEFINE_int32(stage_3_nt, boost::thread::hardware_concurrency(),
    "Total number of parallel threads to use for stage 3");

DEFINE_int32(output_nt, 4,
    "Total number of parallel threads to use for output stage");

DEFINE_int32(output_flag, 1, 
    "Flag to specify output format: "
    "0: BAM (compressed); 1: BAM (uncompressed); 2: SAM");

DEFINE_int32(num_buckets, 1024, 
    "Set output bucket number");

DEFINE_bool(disable_markdup, false,
    "Enable mark duplicate during alignment, default true");

DEFINE_bool(disable_bucketsort, false,
    "Enable bucket sort instead of normal sort, default true");

DEFINE_string(temp_dir, "/tmp", 
    "Temp dir to store intermediate bucket bams");

DEFINE_string(output, "", "Path to output file");

// deprecated
DEFINE_string(output_dir, "",
    "(deprecated) Please use --output instead");

DEFINE_int32(max_batch_records, 40, 
    "(deprecated) Flag to specify how many batch to buffer before sort");

DEFINE_bool(offload, true,
    "(deprecated) Use three compute pipeline stages to enable offloading"
    "workload to accelerators. "
    "If disabled, --use_fpga, --fpga_path will be discard");

DEFINE_string(pac_path, "",
    "(deprecated) File path of the modified reference pac file");

DEFINE_bool(disable_sort, false,
    "disable sorting for output bams");

DEFINE_bool(remove_duplicates, false,
    "remove duplicate reads in bam output");

DEFINE_bool(filter_unmap, false,
    "filter unmapped reads in the output");

DEFINE_bool(merge_bams, true,
    "merge bucket_sort bams");

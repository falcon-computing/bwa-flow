#ifndef CONFIG_H
#define CONFIG_H

#include <gflags/gflags.h>
#include <string>

#include "bwa_wrapper.h"

#define FPGA_RET_PARAM_NUM 5

// GFLAGS parameters
DECLARE_bool(inorder_output);
DECLARE_bool(offload);
DECLARE_bool(M);
DECLARE_bool(sort);
DECLARE_bool(use_fpga);
DECLARE_int32(filter);
DECLARE_int32(chunk_size);
DECLARE_int32(extra_thread);
DECLARE_int32(max_fpga_thread);
DECLARE_int32(max_num_records);
DECLARE_int32(extra_thread);
DECLARE_int32(output_flag);
DECLARE_int32(output_nt);
DECLARE_int32(stage_1_nt);
DECLARE_int32(stage_2_nt);
DECLARE_int32(stage_3_nt);
DECLARE_int32(t);
DECLARE_int32(max_batch_records);
DECLARE_string(fpga_path);
DECLARE_string(pac_path);
DECLARE_string(output_dir);
DECLARE_int32(output_nt);
DECLARE_string(R);

// Global parameters
extern ktp_aux_t* aux;
extern int mpi_rank;
extern int mpi_nprocs;

#endif

#ifndef CONFIG_H
#define CONFIG_H

#include <gflags/gflags.h>
#include "bwa_wrapper.h"

#define FPGA_RET_PARAM_NUM 5

// GFLAGS parameters
DECLARE_bool(inorder_output);
DECLARE_bool(offload);
DECLARE_bool(use_fpga);
DECLARE_int32(chunk_size);
DECLARE_int32(max_fpga_thread);
DECLARE_int32(nt);
DECLARE_int32(stage_1_nt);
DECLARE_int32(stage_2_nt);
DECLARE_int32(stage_3_nt);
DECLARE_string(fpga_path);
DECLARE_string(output_dir);

// Global parameters
extern ktp_aux_t* aux;
#ifdef SCALE_OUT
extern int mpi_rank;
extern int mpi_nprocs;
#endif

#endif

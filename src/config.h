#ifndef CONFIG_H
#define CONFIG_H

#include <gflags/gflags.h>

#define FPGA_RET_PARAM_NUM 5

// GFLAGS parameters
DECLARE_bool(offload);
DECLARE_bool(use_fpga);
DECLARE_string(fpga_path);
DECLARE_int32(chunk_size);
DECLARE_bool(inorder_output);
DECLARE_string(output_dir);
DECLARE_int32(nt);
DECLARE_int32(stage_1_nt);
DECLARE_int32(stage_2_nt);
DECLARE_int32(stage_3_nt);

#endif

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
DECLARE_bool(disable_smem_fpga);
DECLARE_bool(disable_smem_cpu);
DECLARE_bool(disable_sw_fpga);
DECLARE_bool(disable_sw_cpu);
DECLARE_int32(filter);
DECLARE_int32(chunk_size);
DECLARE_int32(max_fpga_thread);
DECLARE_int32(max_num_records);
DECLARE_int32(extra_thread);
DECLARE_int32(output_flag);
DECLARE_int32(stage_1_nt);
DECLARE_int32(stage_2_nt);
DECLARE_int32(stage_3_nt);
DECLARE_int32(t);
DECLARE_string(fpga_path);
DECLARE_string(pac_path);
DECLARE_string(output_dir);
DECLARE_int32(output_nt);
DECLARE_string(R);
DECLARE_int32(max_batch_records);

DECLARE_int32(k);
DECLARE_bool(1);
DECLARE_string(x);
DECLARE_int32(w);
DECLARE_int32(A);
DECLARE_int32(B);
DECLARE_int32(T);
DECLARE_int32(U);
DECLARE_bool(P);
DECLARE_bool(a);
DECLARE_bool(p);
DECLARE_bool(S);
DECLARE_bool(Y);
DECLARE_bool(V);
DECLARE_int32(c);
DECLARE_int32(d);

#if 0
//Overrided by FLAGS_v from glog
DECLARE_bool(v)
#endif

DECLARE_bool(j);
DECLARE_double(r);
DECLARE_double(D);
DECLARE_int32(m);
DECLARE_int32(s);
DECLARE_int32(G); // no description in original bwa codes
DECLARE_int32(N); // no description in original bwa codes
DECLARE_int32(W);
DECLARE_int32(y);
DECLARE_bool(C);
DECLARE_int32(K); // no description in original bwa codes
DECLARE_double(X); // no description in original bwa codes

#ifdef USE_HTSLIB
DECLARE_int32(o);
#endif

DECLARE_string(h);
DECLARE_int32(Q);
DECLARE_string(O);
DECLARE_string(E);
DECLARE_string(L);
DECLARE_string(H);
DECLARE_string(I);

//mark_dup
DECLARE_bool(disable_markdup);

//bucket_sort
DECLARE_int32(num_buckets);

//if use bucket sort
DECLARE_bool(disable_bucketsort);

//for sort_merge flow
DECLARE_string(temp_dir);
DECLARE_string(output);
DECLARE_bool(merge_bams);

DECLARE_bool(disable_sort);
DECLARE_bool(remove_duplicates);
DECLARE_bool(filter_unmap);
// Global parameters
extern ktp_aux_t* aux;
extern int mpi_rank;
extern int mpi_nprocs;

#endif

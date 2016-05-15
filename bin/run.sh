#!/bin/bash

option="--logtostderr=1 \
        --v=2 \
        --offload \
        --use_fpga \
        --fpga_path=/curr/yaoh/kernels/smithwaterman_80pe_new.xclbin \
        --stage_1_nt=7 \
        --stage_2_nt=4 \
        --stage_3_nt=2 \
        --inorder_output"

LD_LIBRARY_PATH=/curr/diwu/tools/gflags/build/lib:$LD_LIBRARY_PATH \
./bwa mem \
    $option \
    /space/scratch/genome/ref/human_g1k_v37.fasta $@


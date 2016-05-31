#!/bin/bash
source ../setup.sh

option="--alsologtostderr=1 \
        --log_dir=. \
        --v=2 \
        --offload \
        --inorder_output"

./bwa mem \
    $option \
    /space/scratch/genome/ref/factor4/human_g1k_v37.fasta $@

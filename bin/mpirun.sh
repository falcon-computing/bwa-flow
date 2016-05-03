#!/bin/bash

MPI_DIR=/curr/diwu/tools/openmpi-1.10.2/build/install

$MPI_DIR/bin/mpirun \
    -np 4 \
    --mca btl_tcp_if_include eth1 \
    --map-by NUMA:PE=6 \
    --bind-to core \
    --host falcon2,falcon3,falcon2,falcon3 \
    ./bwa mem \
    /space/scratch/genome/ref/human_g1k_v37.fasta \
    $@

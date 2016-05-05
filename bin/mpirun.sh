#!/bin/bash

MPI_DIR=/curr/diwu/tools/openmpi-1.10.2/build/install

num_proc=2
host_list=falcon2,falcon2
#host_list=falcon2,falcon3,falcon4,falcon5,falcon2,falcon3,falcon4,falcon5

$MPI_DIR/bin/mpirun \
    -np $num_proc \
    --map-by NUMA:PE=6 \
    --bind-to core \
    --host $host_list \
    --mca btl_tcp_if_include eth1 \
    --mca btl_tcp_sndbuf 1232896 \
    --mca btl_tcp_rcvbuf 1232896 \
    --mca btl_tcp_links 4 \
    ./bwa mem \
    /space/scratch/genome/ref/human_g1k_v37.fasta \
    $@


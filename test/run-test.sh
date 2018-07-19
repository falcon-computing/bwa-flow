#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
source $CURR_DIR/global.bash

$CURR_DIR/data/get-data.sh

$bats_bin $CURR_DIR/env.bats

samples=( $(ls -1 $data_dir/*.fastq.gz | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/') )
for sample in "${samples[@]}"; do
  sample_name=$(basename $sample)
  sample_name=$(echo ${sample_name/_sampled_/})
  n_char=$(( ${#sample_name} + 20 ))
  printf '=%.0s' $(seq 1 $n_char)
  echo ""
  echo "| Testing sample: $sample_name |"
  printf '=%.0s' $(seq 1 $n_char)
  echo ""

  # run unit tests
  GLOG_log_dir=$temp_dir/ \
  LD_LIBRARY_PATH=$OPENMPI_DIR/lib:$LD_LIBRARY_PATH \
  $test_bin mem $ref_genome \
  ${sample}1.fastq.gz \
  ${sample}2.fastq.gz

  # run BATS tests
  TEST_SAMPLE_PREFIX=$sample $bats_bin $CURR_DIR/test.bats
done

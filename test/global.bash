#!/bin/bash
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

local_dir=/local
data_dir=$SCRIPT_DIR/data
temp_dir=$local_dir/${USER}-bwa-test
bats_bin=$SCRIPT_DIR/bats/bin/bats
test_bin=$SCRIPT_DIR/bin/bwa-test
bwa_bin=$SCRIPT_DIR/../bin/bwa
sambamba=$data_dir/sambamba

source $SCRIPT_DIR/test-config.sh

ref_genome=$REF_GENOME

run_cpu() {
  local out=$3;
  local opts=$4;
  $bwa_bin mem \
    --log_dir=$temp_dir \
    --v=2 \
    --output_dir=$out \
    $opts \
    $ref_genome $1 $2;
}

run_fpga() {
  local out=$3;
  local opts=$4;
  $bwa_bin mem \
    --log_dir=$temp_dir \
    --v=2 \
    --chunk_size=2000 \
    --use_fpga \
    --fpga_path=$data_dir/bit.awsxclbin \
    --max_fpga_thread=1 \
    --extra_thread=1 \
    --output_dir=$out \
    $opts \
    $ref_genome $1 $2 ;
}

get_flagstat() {
  local input=$1;
  local output=$2;
  if [ ! -d $input ]; then
    echo "cannot find $input"
    exit 1
  fi;
  local bam=${input}_merged.bam;
  local ret=0;
  if [ "$(ls -1 $input/part-* | wc -l)" -eq 1 ]; then
    cp $input/part-000000 $bam
    ret=$?
  else
    $sambamba merge $bam $(ls $input/part-*)
    ret=$?
  fi;
  if [ $ret -ne 0 ]; then
    echo "failed to get $bam"
    exit 2
  fi;
  $sambamba flagstat $bam > $output;
}

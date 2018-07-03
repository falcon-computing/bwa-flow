#!/bin/bash

get_abs_parentdir() {
  # $1 : relative filename
  filename=$1
  parentdir=$(dirname "${filename}")
  echo "$(cd "${parentdir}" && pwd)"
}

export test_bin=$1
if [[ ! -f $test_bin ]]; then
  echo "$#"
  echo "test file \"$test_bin\" not found!!!"
  exit 1
fi

test_bin_abs_parentdir=$(get_abs_parentdir $test_bin)
bitstream_abs_parentdir=$(get_abs_parentdir $test_bin_abs_parentdir)"/data"
echo "mkdir -p $bitstream_abs_parentdir"
mkdir -p $bitstream_abs_parentdir
echo "aws s3 cp s3://fcs-genome-data/data-suite/bwa-flow/bit.awsxclbin $bitstream_abs_parentdir"
aws s3 cp s3://fcs-genome-data/data-suite/bwa-flow/bit.awsxclbin $bitstream_abs_parentdir


export ref_genome=/genome/ref/v37/human_g1k_v37.fasta
export fastq1=/genome/fastq/sampled/A15_sample_1.fastq.gz
export fastq2=/genome/fastq/sampled/A15_sample_2.fastq.gz

echo "GLOG_v=3 GLOG_logtostderr=1 $test_bin mem $ref_genome $fastq1 $fastq2"
GLOG_v=3 GLOG_logtostderr=1 $test_bin mem $ref_genome $fastq1 $fastq2

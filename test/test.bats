#!/usr/bin/env bats

load global

sample=$TEST_SAMPLE_PREFIX
sample_name=$(basename $sample)
sample_name=$(echo ${sample_name/_sampled_/})

@test "check input sample" {
 [ ! -z "$sample" ]
 [ -f ${sample}1.fastq.gz ]
 [ -f ${sample}2.fastq.gz ]
}

@test "run cpu version" {
  run run_cpu \
    ${sample}1.fastq.gz \
    ${sample}2.fastq.gz \
    $temp_dir/${sample_name}_cpu &> $temp_dir/${sample_name}_cpu.log

  [ $status -eq 0 ]
}

@test "run flagstat on cpu input" {
  run get_flagstat \
    $temp_dir/${sample_name}_cpu \
    $temp_dir/${sample_name}_cpu.flagstat

  [ -f $temp_dir/${sample_name}_cpu.flagstat ] 
}

@test "run fpga version" {
  run run_fpga \
    ${sample}1.fastq.gz \
    ${sample}2.fastq.gz \
    $temp_dir/${sample_name}_fpga &> $temp_dir/${sample_name}_fpga.log

  [ $status -eq 0 ]
}

@test "run flagstat on fpga input" {
  run get_flagstat \
    $temp_dir/${sample_name}_fpga \
    $temp_dir/${sample_name}_fpga.flagstat
  
  [ -f $temp_dir/${sample_name}_fpga.flagstat ]
}

@test "check flagstat" {
  run diff \
    $temp_dir/${sample_name}_cpu.flagstat \
    $temp_dir/${sample_name}_fpga.flagstat

  [ $status -eq 0 ]
}


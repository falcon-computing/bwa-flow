#!/bin/bash

file=$1

check_input() {
  exp=$1;
  col=$2;
  res=`grep "$exp" $file | head`;
  if [[ "$res" != "" ]]; then
    total_time=`grep "$exp" $file | sum.sh $col | sed -e 's/[eE]+*/\\*10\\^/'`;
    total_num=`grep "$exp" $file | wc -l`;
    echo `bc -l <<< "scale=3; $total_time / $total_num"`;
  fi;
}

echo "Read seqs: "`check_input "Read" 1`
echo "Sending seq: "`check_input "Sending seqs batch" 1`
echo "seq2chain: "`check_input "Produced a chain" 1`
echo "prepare ref: "`check_input "prepareChainRef" 1`
echo "SW kernel: "`check_input "SW-FPGA" 1`
echo "Wait for kernel: "`check_input "Wait for FPGA" 1`
echo "FPGA output: "`check_input "FPGA output" 1`
echo "nextTask: "`check_input "nextTask" 4`
echo "packData: "`check_input "packData" 1`
echo "SWRead::Successful: "`check_input "SWRead::Successful" 1`
echo "SWRead::Pending: "`check_input "SWRead::Pending" 1`
echo "SWRead::Finish: "`check_input "SWRead::Finish" 1`
echo "Batch: "`check_input "Batch takes" 1`
echo "FPGA utilization: "`check_input "FPGA utilization" 1`
#echo "Batch: "`check_input "Extension task" 1`
echo "Read: "`check_input "Finished read" 1`
echo "chain2reg: "`check_input "Produced a region" 1`
echo "seedcoverage: "`check_input "Seed coverage" 1`
echo "reg2sam: "`check_input "Produced a sam batch" 1`
echo "seq2sam: "`check_input "Compute a batch" 1`
echo "Sending sam: "`check_input "Sending sam batch" 1`
echo "Output sam: "`check_input "Written" 1`


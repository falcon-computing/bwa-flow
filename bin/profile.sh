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
echo "seq2chain: "`check_input "Finished SeqToChains()" 1`
echo "SW kernel: "`check_input "kernel" 1`
echo "Prepare FPGA data: "`check_input "nextTask" 13`
echo "Process FPGA output: "`check_input "Process output" 1`
echo "Batch: "`check_input "Batch takes" 1`
echo "Wait for input for FPGA: "`check_input "Wait for input for FPGA takes" 1`
echo "chain2reg on CPU: "`check_input "Finished ChainsToRegions() on CPU" 1`
echo "chain2reg on FPGA: "`check_input "Finished ChainsToRegions() on FPGA" 1`
echo "chain2reg average: "`check_input "Finished ChainsToRegions() " 1`
echo "reg2sam: "`check_input "Finished RegionsToSam()" 1`
echo "seq2sam: "`check_input "Compute a batch" 1`
echo "Output sam: "`check_input "Written" 1`
echo "FPGA util: "`check_input "utilization" 1`"%"
echo "Real time: "`check_input "Real time:" 5`" sec"

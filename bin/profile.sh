#!/bin/bash

file=$1

total_time() {
  exp=$1;
  col=$2;
  usec=`grep -i "$exp" $file | sum.sh $col`;
  echo `bc <<< "scale=6; $usec / 1000000"`;
}

total_num() {
  exp=$1;
  grep -i "$exp" $file | wc -l;
}

check_input() {
  exp=$1;
  col=$2;
  res=`grep "$exp" $file | head`;
  if [[ "$res" != "" ]]; then
    ttime=`total_time "$exp" $col`
    tnum=`total_num "$exp"`
    echo `bc -l <<< "scale=3; $ttime / $tnum"`
  fi;
}

#echo "total: "`total_time "Finished" 1`
#echo "seq2chain: "`total_time "Finished SeqToChains()" 1`
#echo "chain2reg: "`total_time "Finished ChainsToRegions()" 1`
#echo "chain2reg on CPU: "`total_time "Finished ChainsToRegions() on CPU" 1`
#echo "chain2reg on CPU: "`total_num "Finished ChainsToRegions() on CPU" 1`
#echo "chain2reg on FPGA: "`total_time "Finished ChainsToRegions() on FPGA" 1`
#echo "chain2reg on FPGA: "`total_num "Finished ChainsToRegions() on FPGA" 1`
#echo "reg2sam: "`total_time "Finished RegionsToSam()" 1`
#echo "FPGA util: "`check_input "utilization" 1`"%"


col_0=`total_time "Read" 1`
col_1=`total_time "Finished SeqsToChains()" 1`
col_2=`total_time "Finished ChainsToRegions()" 1`
col_3=`total_time "Finished RegionsToSam()" 1`
col_4=`total_time "Finished" 1`
col_5=`total_time "Sort" 1`
col_6=`total_time "Written" 1`
col_7=`cat $file | grep "Started\|Finished" | ./perf.py | head -n 1 | awk '{print $NF}'`
col_8=`total_num "Read"`
#col_6=`total_num "Finished ChainsToRegions() on FPGA" 1`
#col_7=`total_num "Finished ChainsToRegions() on CPU" 1`
#col_8=`total_time "Finished ChainsToRegions() on FPGA" 1`
#col_9=`total_time "Finished ChainsToRegions() on CPU" 1`

printf "%f, %f, %f, %f, %f, %f, %f, %f, %d\n" \
        $col_0 $col_1 $col_2 $col_3 $col_4 $col_5 $col_6 $col_7 $col_8

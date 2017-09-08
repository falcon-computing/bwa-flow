#!/bin/bash

pid=$1
output=$(ps -o pid,pcpu,pmem $pid | sed -n 2p)

while [ ! -z "$output" ]; do
  output=$(ps -o pid,pcpu,pmem $pid | sed -n 2p)
  pcpu=$(echo $output | awk '{print $2}')
  pmem=$(echo $output | awk '{print $3}')
  printf "%s, %s, %s\n" "$(date +"%F %T")" "$pcpu" "$pmem"
  sleep 10
done

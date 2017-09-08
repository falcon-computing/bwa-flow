#!/bin/bash
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

profile=default
bwa=../bin/bwa
bwa_options=

if [ $# -gt 0 ]; then
  profile=$1
  if [ -f ${SCRIPT_DIR}/profile_${profile}.sh ]; then
    source ${SCRIPT_DIR}/profile_${profile}.sh
    echo "Use profile: $profile"
  else
    echo "Cannot find profile $profile, use default settings"
    profile=default
  fi
fi

run_id=${profile}"-"$(date +%F-%T)
log_dir=${run_id}
mkdir -p $log_dir

function get_sample_name() {
  local dir=$1;
  ls $dir/*.fastq.gz | xargs -i basename {} | sed -e 's/\.fastq\.gz$//g' | sed -e 's/_[12]//g' | sort -u;
}

for sample in $(get_sample_name $input_dir); do
  $bwa mem \
    --v=1 \
    --logtostderr=1 \
    --output_dir=$tmp_dir/$sample \
    $bwa_options \
    $ref_genome \
    $(ls $input_dir/${sample}*) 1>/dev/null 2> $log_dir/${sample}.out &

  pid=$!
  ./monitor.sh $pid > $log_dir/${sample}.perf &

  wait $pid

  if [ $? -ne 0 ]; then
    echo "sample $sample failed"
    echo "$sample, -1" >> perf_${run_id}.log
  else
    e_time=$(cat $log_dir/${sample}.out | grep "Real time" | awk '{print $(NF-5)}')
    echo "finished sample $sample in $e_time sec"
    echo "$sample, $e_time" >> perf_${run_id}.log
  fi

  rm -rf $tmp_dir/$sample
done

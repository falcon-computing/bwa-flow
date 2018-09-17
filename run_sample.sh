#!/bin/bash

############################################################
#   Global setting
############################################################

export GLOG_v=2
export GLOG_logtostderr=1

export TEST_ONE=1
export TEST_MODE=2


############################################################
#   Test cases setting
############################################################

export REF_GENOME=/local/ref/human_g1k_v37.fasta
#export REF_GENOME=/local/ref/human_g1k_v37.fasta
#export REF_GENOME=/genome/ref/human_g1k_v37.fasta
export FASTQ_DIR=/local/fastq/WES
export FASTQ_REX="^\(.*\)-.*_R[12].*.fastq.gz$/\1"
export OUTPUT_DIR=/pool/storage/jyqiu/tests


############################################################
#   Build setting
############################################################

export EXE=$(pwd)/bin/bwa
#export BIT_STREAM="--sw_fpga_path="$(pwd)/mem_collect_intv_sw_top_xilinx_vcu1525_dynamic_5_0_0906.xclbin
export BIT_STREAM="--sw_fpga_path="$(pwd)/mem_collect_intv_sw_top_xilinx_vcu1525_dynamic_5_0_0910.xclbin
#export BIT_STREAM=$(pwd)/mem_collect_intv.xclbin
#export BIT_STREAM="--smem_fpga_path="$(pwd)/mem_collect_intv_xilinx_vcu1525_dynamic_5_0_0827.xclbin


############################################################
#   Running body
############################################################

export file_list_=( $(ls ${FASTQ_DIR}) )
export case_list_=()
for f in ${file_list_[@]};
do
  export c=$(echo $f | sed 's/'${FASTQ_REX}'/g' | sed 's/.fastq.gz//g' | sed 's/.fastq//g')
  if [[ ! " ${case_list_[@]} " =~ " $c " ]];
  then
    case_list_+=( $c )
  fi
done

for c in ${case_list_[@]};
do
  export fastq_files_=()
  for f in ${file_list_[@]};
  do
    if [[ $f =~ ^.*${c}.*$ ]];
    then
      fastq_files_+=( $f )
    fi
  done
  
  if [[ ${#fastq_files_[@]} -eq 1 ]];
  then
     export FASTQ1=${FASTQ_DIR}/${fastq_files_[0]}
     export FASTQ2=
  elif [[ ${#fastq_files_[@]} -eq 2 ]]; 
  then
     export FASTQ1=${FASTQ_DIR}/${fastq_files_[0]}
     export FASTQ2=${FASTQ_DIR}/${fastq_files_[1]}
  else
     echo [$c]
     echo "ERROR: "${fastq_files_[*]}
     continue
  fi

  #HACK
  export c=Garvan_sampled
  export FASTQ1=/pool/storage/fastq/NA12878-Garvan-Vial1_R1.sampled.fastq.gz
  export FASTQ2=/pool/storage/fastq/NA12878-Garvan-Vial1_R2.sampled.fastq.gz
  export TEST_ONE=1
  

  export CASE_NAME=$c
  export DATE_TAG=$(date +%y_%m_%d_%H_%M_%S)
  export TEST_DIR=${OUTPUT_DIR}/${CASE_NAME}_${DATE_TAG}
  mkdir -p ${TEST_DIR}
  export BAM_DIR=${TEST_DIR}/bam
  export LOG_DIR=${TEST_DIR}/log
  mkdir -p ${LOG_DIR}

  #echo ${CASE_NAME}:${FASTQ1}:${FASTQ2}
  if [[ $TEST_MODE -eq 1 || $TEST_MODE -eq 3 ]]; then
    echo "[${CASE_NAME}:CPU]"
    $EXE mem \
      --log_dir=${LOG_DIR} \
      --v=${GLOG_v} \
      --output_dir=${BAM_DIR} \
      ${REF_GENOME} \
      ${FASTQ1} \
      ${FASTQ2} \
    2>&1 | tee ${TEST_DIR}/xilinx_cpu.log
  fi
  
  if [[ $TEST_MODE -eq 2 || $TEST_MODE -eq 3 ]]; then
    echo "[${CASE_NAME}:FPGA]"
    $EXE mem \
      --log_dir=${LOG_DIR} \
      --v=${GLOG_v} \
      --chunk_size=2000 \
      --use_fpga \
      ${BIT_STREAM} \
      --inorder_output \
      --extra_thread=0 \
      --output_dir=${BAM_DIR} \
      ${REF_GENOME} \
      ${FASTQ1} \
      ${FASTQ2} \
    2>&1 | tee ${TEST_DIR}/xilinx_fpga.log
  fi
  
  if [[ $TEST_ONE == 1 ]]; then
    break
  fi

  echo

done

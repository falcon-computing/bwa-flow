#!/bin/bash
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
s3_dir=s3://fcs-data/gatk3-build/deps/$(uname -s)

cd $DIR
aws s3 cp --recursive --include "*.tar.gz" $s3_dir/ .

for tar in `ls *.tar.gz`; do
  tar zxf $tar 
  if [ $? -ne 0 ]; then
    cd -
    exit -1
  fi
  rm $tar
done

echo "ready" > .ready
cd -

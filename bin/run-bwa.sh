#!/bin/bash
BWA_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

source $BWA_DIR/../setup.sh
source /curr/diwu/setenv_sdaccel2015.4.sh

option="--logtostderr=1 \
        --v=1 \
        --offload \
        --nt=32"

$BWA_DIR/bwa $option $@

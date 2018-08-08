#!/bin/bash
# This preparation file helps to yield compatibility with previous config.mk file.

CURDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# basic reference file
REF_GENOME=/local/ref/human_g1k_v37.fasta
# bwa-flow executable
BWA_BIN=${CURDIR}/../build/bwa-flow
# test executable
BWA_TEST_BIN=${CURDIR}/../build/test/test_bwa

# set up procedure
mkdir -p ${CURDIR}/../bin
cp ${BWA_BIN} ${CURDIR}/../bin/bwa
mkdir -p ${CURDIR}/bin
cp ${BWA_TEST_BIN} ${CURDIR}/bin/bwa-test

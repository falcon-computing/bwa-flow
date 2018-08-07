#!/bin/bash
CURR_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)
aws s3 sync s3://fcs-genome-data/data-suite/bwa-flow/ $CURR_DIR/

#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
source $DIR/config.mk

export LD_LIBRARY_PATH=$BOOST_DIR/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GLOG_DIR/lib:$LD_LIBRARY_PATH

BWA_DIR   := ../bwa
KFLOW_DIR := ../kflow
SRC_DIR   := ../src

#OPENMPI_DIR=/curr/diwu/tools/openmpi-1.10.2/build/install
include ../config.mk

CC	:= gcc
PP	:= g++

CFLAGS := -c -fPIC -std=c++0x -g -O2 

INCLUDES := -I.. -I../src \
	    -I$(KFLOW_DIR)/include \
	    -I$(BOOST_DIR)/include \
	    -I$(GLOG_DIR)/include \
	   -I$(GFLAGS_DIR)/include \
	    -I$(GTEST_DIR)/include


LINK    := -L$(BWA_DIR) -lbwa \
	   -L$(KFLOW_DIR)/lib -lkflow \
	   -L$(BOOST_DIR)/lib \
	   	-lboost_system \
	   	-lboost_thread \
	   	-lboost_iostreams \
	   	-lboost_filesystem \
	   	-lboost_regex \
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GFLAGS_DIR)/lib -lgflags \
	   -L$(GTEST_DIR)/build -lgtest \
	   -lpthread -lm -ldl -lz -lrt

ifneq ($(HTSLIB_PATH),)
CFLAGS   := $(CFLAGS) -DUSE_HTSLIB
INCLUDES := $(INCLUDES) -I$(HTSLIB_PATH)
LINK     := $(LINK) -L$(HTSLIB_PATH) -lhts 
endif 

ifneq ($(OPENMPI_DIR),)
INCLUDES := $(INCLUDES) -I$(OPENMPI_DIR)/include
LINK	 := $(LINK) -L$(OPENMPI_DIR)/lib -lmpi_cxx -lmpi
endif

COMPILE := $(CFLAGS) $(INCLUDES)

OBJS   	:= main.o \
	   CommTest.o
DEP_OBJS:= $(SRC_DIR)/Pipeline.o \
	   $(SRC_DIR)/MPIPipeline.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/util.o

REF_GENOME  := /pool/local/ref/human_g1k_v37.fasta
TEST_FASTQ1 := /pool/storage/fastq/gcat_set_041_1.fastq
TEST_FASTQ2 := /pool/storage/fastq/gcat_set_041_2.fastq

DST	:= test_app

.PHONY: all test clean

all: $(DST)

test: 
	GLOG_v=3 \
	GLOG_logtostderr=1 \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	./$(DST) mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)


$(DST): $(OBJS) $(DEP_OBJS)
	$(PP) -o $@ $(OBJS) $(DEP_OBJS) $(LINK) 

$(DEP_OBJS): 
	make -C ../

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -f $(OBJS)
	rm -f $(DST)


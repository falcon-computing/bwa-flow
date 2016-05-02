BWA_DIR   := ../bwa
KFLOW_DIR := ../kflow
SRC_DIR   := ../src

include /curr/diwu/prog/blaze/Makefile.config

CC	:= gcc
PP	:= g++

CFLAGS := -c -fPIC -std=c++0x -g

INCLUDES := -I.. -I../src \
	    -I$(KFLOW_DIR)/include \
	    -I$(BOOST_DIR)/include \
	    -I$(GLOG_DIR)/include \
	    -I$(GTEST_DIR)/include \
	    -I$(XILINX_OPENCL_DIR)/runtime/include/1_2

COMPILE := $(CFLAGS) $(INCLUDES)

LINK    := -L$(BWA_DIR) -lbwa \
	   -L$(KFLOW_DIR)/lib -lkflow \
	   -L$(BOOST_DIR)/lib \
	   	-lboost_system \
	   	-lboost_thread \
	   	-lboost_iostreams \
	   	-lboost_filesystem \
	   	-lboost_regex \
	   -L$(PROTOBUF_DIR)/lib -lprotobuf \
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GTEST_DIR)/build -lgtest \
	   -L$(XILINX_OPENCL_DIR)/runtime/lib/x86_64 -lOpenCL \
	   -lpthread -lm -ldl

OBJS   	:= main.o \
	   CommTest.o \
           $(SRC_DIR)/Extension.o\
           $(SRC_DIR)/FPGAAgent.o \
	   $(SRC_DIR)/Pipeline.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/SWRead.o \
	   $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/util.o

DST	:= test_app

.PHONY: all test clean

all: $(DST)

test: 
	GLOG_v=3 \
	GLOG_logtostderr=1 \
	./$(DST) mem \
		/space/scratch/genome/ref/human_g1k_v37.fasta \
		/space/scratch/genome/fastq/gcat_set_041_1.fastq \
		/space/scratch/genome/fastq/gcat_set_041_2.fastq


$(DST): $(OBJS)
	$(PP) -o $@ $(OBJS) $(LINK) 

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -f $(OBJS)
	rm -f $(DST)


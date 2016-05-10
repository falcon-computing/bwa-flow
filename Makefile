include /curr/diwu/prog/blaze/Makefile.config

BWA_DIR   	:= ./bwa
KFLOW_DIR 	:= ./kflow
SRC_DIR   	:= ./src
MANAGER_DIR 	:= /curr/diwu/prog/blaze/manager
GFLAGS_DIR	:= /curr/diwu/tools/gflags/build
OPENMPI_DIR	:= /curr/diwu/tools/openmpi-1.10.2/build/install

MAKE	:= make
CC	:= gcc
PP	:= g++
MPIPP 	:= $(OPENMPI_DIR)/bin/mpic++

CFLAGS 	:= -g -std=c++0x -fPIC -O3
OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/Pipeline.o \
           $(SRC_DIR)/Extension.o\
           $(SRC_DIR)/FPGAAgent.o \
	   $(SRC_DIR)/main.o \
	   $(SRC_DIR)/SWRead.o \
	   $(SRC_DIR)/util.o

INCLUDES:= -I. -I$(BWA_DIR) \
	   -I$(KFLOW_DIR)/include \
	   -I$(BOOST_DIR)/include \
	   -I$(XILINX_OPENCL_DIR)/runtime/include/1_2 \
	   -I$(PROTOBUF_DIR)/include \
	   -I$(GLOG_DIR)/include \
	   -I$(GFLAGS_DIR)/include
	
LIBS	:= -L$(BWA_DIR) -lbwa \
	   -L$(KFLOW_DIR)/lib -lkflow \
	   -L$(BOOST_DIR)/lib \
	   	-lboost_system \
		-lboost_thread \
		-lboost_iostreams \
		-lboost_filesystem \
		-lboost_regex \
	   -L$(PROTOBUF_DIR)/lib -lprotobuf \
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GFLAGS_DIR)/lib -lgflags \
	   -L$(XILINX_OPENCL_DIR)/runtime/lib/x86_64 -lOpenCL \
	   -lpthread -lm -ldl -lz -lrt

ifneq ($(NDEBUG),)
CFLAGS   := $(CFLAGS) -DNDEBUG
endif

ifneq ($(SCALEOUT),)
CFLAGS   := $(CFLAGS) -DSCALE_OUT
INCLUDES := $(INCLUDES) -I$(OPENMPI_DIR)/include
LIBS	 := $(LIBS) -L$(OPENMPI_DIR)/lib -lmpi_cxx -lmpi
endif

PROG	 := ./bin/bwa

all:	$(PROG)

scaleout:
	$(MAKE) SCALEOUT=1 all

./bin/bwa: $(BWA_DIR)/libbwa.a $(OBJS) 
	$(PP) $(OBJS) -o $@ $(LIBS)

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

./bwa/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f $(OBJS)
	rm -f $(PROG)  

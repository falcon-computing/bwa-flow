MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR  := $(dir $(MKFILE_PATH))
include $(MKFILE_DIR)/config.mk

BIN_DIR   	:= $(MKFILE_DIR)/bin
BWA_DIR   	:= $(MKFILE_DIR)/bwa
KFLOW_DIR 	:= $(MKFILE_DIR)/kflow
SRC_DIR   	:= $(MKFILE_DIR)/src
TEST_DIR	:= $(MKFILE_DIR)/test

CFLAGS 	:= -std=c++0x -fPIC -O3 

OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/Pipeline.o \
	   $(SRC_DIR)/util.o

STDOBJS := $(SRC_DIR)/main.o 

MPIOBJS := $(SRC_DIR)/MPIPipeline.o \
	   $(SRC_DIR)/mpi_main.o

TESTOBJS:= $(TEST_DIR)/main.o \
	   $(TEST_DIR)/PipelineTests.o \
	   $(TEST_DIR)/UtilTests.o

TEST_DEPOBJS := $(SRC_DIR)/Pipeline.o \
	   	$(SRC_DIR)/MPIPipeline.o \
	   	$(SRC_DIR)/preprocess.o \
	   	$(SRC_DIR)/wrappered_mem.o \
	   	$(SRC_DIR)/util.o


INCLUDES:= -I$(MKFILE_DIR) \
	   -I$(SRC_DIR) \
	   -I$(KFLOW_DIR)/include \
	   -I$(BOOST_DIR)/include \
	   -I$(GLOG_DIR)/include \
	   -I$(GFLAGS_DIR)/include \
	   -I$(GTEST_DIR)/include
	
LIBS	:= -L$(BWA_DIR) -lbwa \
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

GIT_VERSION := $(shell git describe --abbrev=5 --dirty --always --tags)
CFLAGS	:= $(CFLAGS) -DVERSION=\"$(GIT_VERSION)\"

PROG	 := $(BIN_DIR)/bwa
MPIPROG	 := $(BIN_DIR)/bwa-mpi
TESTPROG := $(TEST_DIR)/bwa-test


ifneq ($(DEBUG),)
CFLAGS   := $(CFLAGS) -g
else
CFLAGS   := $(CFLAGS) -DNDEBUG
endif

ifneq ($(HTSLIB_PATH),)
CFLAGS   := $(CFLAGS) -DUSE_HTSLIB
INCLUDES := $(INCLUDES) -I$(HTSLIB_PATH)
LIBS     := $(LIBS) -L$(HTSLIB_PATH) -lhts 
endif 

ifneq ($(BUILD_FPGA),)
CFLAGS 	 := $(CFLAGS) -DBUILD_FPGA
OBJS	 := $(OBJS) \
	    $(SRC_DIR)/FPGAPipeline.o \
            $(SRC_DIR)/FPGAAgent.o 
INCLUDES := $(INCLUDES) \
	    -I$(XILINX_OPENCL_DIR)/runtime/include/1_2 
LIBS	 := $(LIBS) \
	    -L$(XILINX_OPENCL_DIR)/runtime/lib/x86_64 -lOpenCL
endif

ifneq ($(OPENMPI_DIR),)
INCLUDES := $(INCLUDES) -I$(OPENMPI_DIR)/include
MPILIBS	 := -L$(OPENMPI_DIR)/lib -lmpi_cxx -lmpi
MPIPROG	 := ./bin/bwa-mpi
endif

# check FLMDIR
ifneq ($(FLMDIR),)
# add support for flex license manage
FLMLIB 		:= -llmgr_trl -lcrvs -lsb -lnoact -llmgr_dongle_stub

CFLAGS   	:= $(CFLAGS) -DNDEBUG -DUSELICENSE
INCLUDES 	:= $(INCLUDES) -I$(FLMDIR)
LIBS		:= $(LIBS) -L$(FLMDIR) $(FLMLIB) 
LMDEPS 	 	:= $(FLMDIR)/license.o \
		   $(FLMDIR)/lm_new.o
endif 

all:	$(PROG) $(TESTPROG)

scaleout: $(MPIPROG)

test:	$(TESTPROG)

runtest: 
	GLOG_v=3 \
	GLOG_alsologtostderr=1 \
	GLOG_log_dir=$(TEST_DIR) \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	$(TESTPROG) mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)

$(PROG): $(BWA_DIR)/libbwa.a $(OBJS) $(STDOBJS) $(LMDEPS)
	$(PP) $(OBJS) $(STDOBJS) $(LMDEPS) -o $@ $(LIBS)

$(MPIPROG): $(BWA_DIR)/libbwa.a $(MPIOBJS) $(OBJS) $(LMDEPS)
	$(PP) $(OBJS) $(MPIOBJS) $(LMDEPS) -o $@ $(MPILIBS) $(LIBS)

$(TESTPROG): $(TESTOBJS) $(TEST_DEPOBJS)
	$(PP) $(TESTOBJS) $(TEST_DEPOBJS) -o $@ $(MPILIBS) $(LIBS) 

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(BWA_DIR)/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/*.o
	rm -f $(PROG) $(MPIPROG) $(TESTPROG)

.PHONY: all scaleout test runtest clean

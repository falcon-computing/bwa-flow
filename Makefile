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

MPIOBJS := $(SRC_DIR)/MPIChannel.o \
					 $(SRC_DIR)/mpi_main.o

TESTOBJS:= $(TEST_DIR)/main.o \
	   $(TEST_DIR)/PipelineTests.o \
	   $(TEST_DIR)/UtilTests.o

TEST_DEPOBJS := $(SRC_DIR)/Pipeline.o \
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
	   -L$(GTEST_DIR)/lib -lgtest \
	   -lpthread -lm -ldl -lz -lrt

PROG	 := $(BIN_DIR)/bwa
TESTPROG := $(TEST_DIR)/bwa-test

DEPS	 := ./deps/.ready

ifneq ($(RELEASE),)
CFLAGS   := $(CFLAGS) -DNDEBUG
TESTOPT  := 
else
CFLAGS   := $(CFLAGS) -g
TESTOPT  := GLOG_v=3 \
	    GLOG_alsologtostderr=1
endif

ifneq ($(HTSLIB_PATH),)
CFLAGS   := $(CFLAGS) -DUSE_HTSLIB
INCLUDES := $(INCLUDES) -I$(HTSLIB_PATH)
LIBS     := $(LIBS) -L$(HTSLIB_PATH) -lhts 
endif 

GIT_VERSION := $(shell git describe --tags | sed 's/\(.*\)-.*/\1/')

ifneq ($(BUILD_FPGA),)
CFLAGS 	 := $(CFLAGS) -DBUILD_FPGA
OBJS	 := $(OBJS) \
	    $(SRC_DIR)/FPGAPipeline.o \
	    $(SRC_DIR)/SWTask.o

ifneq ($(ALTERAOCLSDKROOT),)
CFLAGS 	 := $(CFLAGS) -DINTEL_FPGA
OBJS	 := $(OBJS) \
	    $(SRC_DIR)/IntelAgent.o

INCLUDES := $(INCLUDES) \
	    $(shell aocl compile-config )
LIBS	 := $(LIBS) \
	    $(shell aocl link-config )
GIT_VERSION := $(GIT_VERSION)-intel
else 

ifeq ($(XILINX_SDX),)
XILINX_SDX := $(XILINX_OPENCL)
endif

CFLAGS 	 := $(CFLAGS) -DXILINX_FPGA
INCLUDES := $(INCLUDES) \
	-I$(XILINX_SDX)/runtime/include/1_2 
LIBS	 := $(LIBS) \
	-L$(XILINX_SDX)/runtime/lib/x86_64 -lxilinxopencl

OBJS	 := $(OBJS) \
	    $(SRC_DIR)/XCLAgent.o

GIT_VERSION := $(GIT_VERSION)-xlnx
endif
endif

ifeq ($(RELEASE),)
GIT_VERSION := $(GIT_VERSION)-dev
endif

CFLAGS	:= $(CFLAGS) -DVERSION=\"$(GIT_VERSION)\"

ifneq ($(OPENMPI_DIR),)
CFLAGS 	 := $(CFLAGS) -DUSE_MPI
INCLUDES := $(INCLUDES) -I$(OPENMPI_DIR)/include
MPILIBS	 := -L$(OPENMPI_DIR)/lib -lmpi_cxx -lmpi
MPIPROG	 := ./bin/bwa-mpi

TESTOBJS:= $(TESTOBJS) \
	   $(TEST_DIR)/ChannelTests.o

TEST_DEPOBJS := $(TEST_DEPOBJS) \
	   	$(SRC_DIR)/MPIChannel.o
endif

# check FLMDIR
ifneq ($(RELEASE),)
ifneq ($(FLMDIR),)
# add support for flex license manage
FLMLIB 		:= -llmgr_pic_trl -lcrvs -lsb -lnoact -llmgr_dongle_stub_pic

CFLAGS   	:= $(CFLAGS) -DUSELICENSE
INCLUDES 	:= $(INCLUDES) -I$(FLMDIR)/include
LIBS		:= $(LIBS) -L$(FLMDIR)/lib $(FLMLIB) 
LMDEPS 	 	:= $(FLMDIR)/lib/license.o \
		   $(FLMDIR)/lib/lm_new.o
endif 
endif

all:	$(PROG) $(TESTPROG) $(MPIPROG)

dist:   $(PROG)
	aws s3 cp bin/bwa s3://fcs-genome-build/bwa/bwa-$(GIT_VERSION)

scaleout: $(MPIPROG)

test:	$(TESTPROG)

runtest: 
	$(TESTOPT) \
	GLOG_log_dir=$(TEST_DIR) \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	$(TESTPROG)  \
	mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)

runmpitest: 
	$(TESTOPT) \
	GLOG_log_dir=$(TEST_DIR) \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	$(OPENMPI_DIR)/bin/mpirun -np 4 \
	--mca orte_base_help_aggregate 0 \
	$(TESTPROG) --gtest_filter=ChannelTests.* \
	mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)

$(DEPS): ./deps/get-all.sh
	./deps/get-all.sh

$(PROG): $(BWA_DIR)/libbwa.a $(OBJS) $(STDOBJS) $(LMDEPS)
	$(PP) $(OBJS) $(STDOBJS) $(LMDEPS) -o $@ $(LIBS)

$(MPIPROG): $(BWA_DIR)/libbwa.a $(MPIOBJS) $(OBJS) $(LMDEPS)
	$(PP) $(OBJS) $(MPIOBJS) $(LMDEPS) -o $@ $(MPILIBS) $(LIBS)

$(TESTPROG): $(TESTOBJS) $(TEST_DEPOBJS)
	$(PP) $(TESTOBJS) $(TEST_DEPOBJS) -o $@ $(MPILIBS) $(LIBS) 

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp $(SRC_DIR)/*.h
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(TEST_DIR)/%.o: $(TEST_DIR)/%.cpp $(SRC_DIR)/*.h
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(BWA_DIR)/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/*.o
	rm -f $(PROG) $(MPIPROG) $(TESTPROG)

.PHONY: all scaleout test runtest clean

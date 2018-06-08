MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
MKFILE_DIR  := $(dir $(MKFILE_PATH))

GLOG_DIR    := ./deps/glog-0.3.3
GFLAGS_DIR  := ./deps/gflags
GTEST_DIR   := ./deps/googletest
FLMDIR      := ./deps/falcon-lic

include $(MKFILE_DIR)/config.mk

BIN_DIR   	:= $(MKFILE_DIR)/bin
BWA_DIR   	:= $(MKFILE_DIR)/bwa
KFLOW_DIR 	:= $(MKFILE_DIR)/kflow
SRC_DIR   	:= $(MKFILE_DIR)/src
TEST_DIR	:= $(MKFILE_DIR)/test

CFLAGS 	:= -std=c++0x -fPIC

OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/config.o \
	   $(SRC_DIR)/Pipeline.o \
	   $(SRC_DIR)/util.o

STDOBJS := $(SRC_DIR)/main.o 

MPIOBJS := $(SRC_DIR)/MPIChannel.o \
	   $(SRC_DIR)/mpi_main.o

TESTOBJS:= $(TEST_DIR)/src/main.o \
	   $(TEST_DIR)/src/PipelineTests.o \
	   $(TEST_DIR)/src/UtilTests.o

TEST_DEPOBJS := $(SRC_DIR)/Pipeline.o \
	   	$(SRC_DIR)/config.o \
	   	$(SRC_DIR)/preprocess.o \
	   	$(SRC_DIR)/wrappered_mem.o \
	   	$(SRC_DIR)/util.o


INCLUDES:= -I$(MKFILE_DIR) \
	   -I$(SRC_DIR) \
	   -I$(KFLOW_DIR)/include \
	   -I/curr/diwu/prog/blaze/build/install/include \
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
	   -L/curr/diwu/prog/blaze/build/install/lib -lblaze -lblaze_message \
	   -lprotobuf \
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GFLAGS_DIR)/lib -lgflags \
	   -L$(GTEST_DIR)/lib -lgtest \
	   -lpthread -lm -ldl -lz -lrt

PROG	 := $(BIN_DIR)/bwa
TESTPROG := $(TEST_DIR)/bin/bwa-test

DEPS	 := ./deps/.ready

ifneq ($(DEBUG),)
CFLAGS   := $(CFLAGS) -g
TESTOPT  := GLOG_v=3 \
	    GLOG_alsologtostderr=1
else
CFLAGS   := $(CFLAGS) -O3
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

TESTOBJS := $(TESTOBJS) \
	    $(TEST_DIR)/src/FPGATests.o

TEST_DEPOBJS := $(TEST_DEPOBJS) \
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

TEST_DEPOBJS := $(TEST_DEPOBJS) \
	    $(SRC_DIR)/XCLAgent.o

GIT_VERSION := $(GIT_VERSION)-xlnx
endif
endif

ifeq ($(RELEASE),)
GIT_VERSION := $(GIT_VERSION)-dev
endif

ifneq ($(OPENMPI_DIR),)
CFLAGS 	 := $(CFLAGS) -DUSE_MPI
INCLUDES := $(INCLUDES) -I$(OPENMPI_DIR)/include
MPILIBS	 := -L$(OPENMPI_DIR)/lib -lmpi_cxx -lmpi
MPIPROG	 := ./bin/bwa-mpi

TESTOBJS:= $(TESTOBJS) \
	   $(TEST_DIR)/src/ChannelTests.o

TEST_DEPOBJS := $(TEST_DEPOBJS) \
	   	$(SRC_DIR)/MPIChannel.o
endif

ifneq ($(RELEASE),)
# add support for falcon license client
CFLAGS   	:= $(CFLAGS) -DNDEBUG -DUSELICENSE
INCLUDES 	:= $(INCLUDES) -I$(FLMDIR)/include
LIBS 	 	:= -L$(FLMDIR)/lib -lfalcon_license \
		   $(LIBS) 
ifneq ($(DEPLOYMENT),) # config license client for a cloud
CFLAGS       := $(CFLAGS) -DDEPLOY_$(DEPLOYMENT)
GIT_VERSION  := $(GIT_VERSION)-$(DEPLOYMENT)
endif
endif

CFLAGS	:= $(CFLAGS) -DVERSION=\"$(GIT_VERSION)\"

all:	$(PROG) $(TESTPROG) $(MPIPROG)

dist:   $(PROG)
	aws s3 cp bin/bwa s3://fcs-genome-build/bwa/bwa-$(GIT_VERSION)

scaleout: $(MPIPROG)

test:	$(TESTPROG)

runtest: 
	$(TESTOPT) \
	GLOG_log_dir=$(TEST_DIR)/bin \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	$(TESTPROG)  \
	mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)

runmpitest: 
	$(TESTOPT) \
	GLOG_log_dir=$(TEST_DIR)/bin \
	LD_LIBRARY_PATH=$(OPENMPI_DIR)/lib:$(LD_LIBRARY_PATH) \
	$(OPENMPI_DIR)/bin/mpirun -np 4 \
	--mca orte_base_help_aggregate 0 \
	$(TESTPROG) --gtest_filter=ChannelTests.* \
	mem $(REF_GENOME) $(TEST_FASTQ1) $(TEST_FASTQ2)

$(DEPS): ./deps/get-all.sh
	./deps/get-all.sh

$(PROG): $(BWA_DIR)/libbwa.a $(OBJS) $(STDOBJS)
	$(PP) $(OBJS) $(STDOBJS) -o $@ $(LIBS)

$(MPIPROG): $(BWA_DIR)/libbwa.a $(MPIOBJS) $(OBJS)
	$(PP) $(OBJS) $(MPIOBJS) -o $@ $(MPILIBS) $(LIBS)

$(TESTPROG): $(TESTOBJS) $(TEST_DEPOBJS)
	$(PP) $(TESTOBJS) $(TEST_DEPOBJS) -o $@ $(MPILIBS) $(LIBS) 

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp $(SRC_DIR)/*.h
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(TEST_DIR)/src/%.o: $(TEST_DIR)/src/%.cpp $(SRC_DIR)/*.h
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(BWA_DIR)/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(TEST_DIR)/src/*.o
	rm -f $(PROG) $(MPIPROG) $(TESTPROG)

.PHONY: all scaleout test runtest clean

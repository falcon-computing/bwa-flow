include config.mk

BWA_DIR   	:= ./bwa
KFLOW_DIR 	:= ./kflow
SRC_DIR   	:= ./src

CFLAGS 	:= -std=c++0x -fPIC -O3 

OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/Pipeline.o \
	   $(SRC_DIR)/util.o

STDOBJS := $(SRC_DIR)/main.o 

MPIOBJS := $(SRC_DIR)/MPIPipeline.o \
	   $(SRC_DIR)/mpi_main.o

INCLUDES:= -I. -I$(BWA_DIR) \
	   -I$(KFLOW_DIR)/include \
	   -I$(BOOST_DIR)/include \
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
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GFLAGS_DIR)/lib -lgflags \
	   -lpthread -lm -ldl -lz -lrt

GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
CFLAGS	:= $(CFLAGS) -DVERSION=\"$(GIT_VERSION)\"

PROG	 := ./bin/bwa

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
INCLUDES := $(INCLUDES) -I$(OPENMPI_INC)
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



all:	$(PROG) $(MPIPROG)

scaleout: $(MPIPROG)

./bin/bwa-mpi: $(BWA_DIR)/libbwa.a $(MPIOBJS) $(OBJS) $(LMDEPS)
	$(PP) $(OBJS) $(MPIOBJS) $(LMDEPS) -o $@ $(MPILIBS) $(LIBS)

./bin/bwa: $(BWA_DIR)/libbwa.a $(OBJS) $(STDOBJS) $(LMDEPS)
	$(PP) $(OBJS) $(STDOBJS) $(LMDEPS) -o $@ $(LIBS)

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

./bwa/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f $(OBJS) 
	rm -f $(STDOBJS)
	rm -f $(MPIOBJS)
	rm -f $(PROG) $(MPIPROG)

.PHONY: all scaleout clean

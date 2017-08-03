BLAZE_HOME=/curr/diwu/prog/blaze
ifeq ("$(wildcard $(BLAZE_HOME)/Makefile.config)","")
$(error BLAZE_HOME not set properly in environment)
endif

MANAGER_DIR=$(BLAZE_HOME)/manager
PLATFORM_DIR=$(BLAZE_HOME)/platforms/xlnx_opencl
include $(BLAZE_HOME)/Makefile.config

PP=$(XILINX_OPENCL_DIR)/lnx64/tools/gcc/bin/g++
CC=$(XILINX_OPENCL_DIR)/lnx64/tools/gcc/bin/gcc

CFLAGS= -shared -fPIC -std=c++0x -m64 -g 
COMPILE= $(CFLAGS) \
	 -I$(MANAGER_DIR)/include \
	 -I$(PLATFORM_DIR)/include \
	 -I$(BOOST_DIR)/include \
	 -I$(PROTOBUF_DIR)/include \
	 -I$(GLOG_DIR)/include \
	 -I$(GFLAGS_DIR)/include \
	 -I$(XILINX_OPENCL_DIR)/runtime/include/1_2

LINK= -L$(XILINX_OPENCL_DIR)/runtime/lib/x86_64 -lxilinxopencl

ACC = SmithWaterman

OBJS= $(addsuffix .so, $(ACC))

all: $(OBJS)

%.so: %.cpp Makefile
	$(PP) $(COMPILE) $(LINK) $< -o $@

clean:
	rm -rf $(OBJS)

ROOT_DIR:= .
include $(ROOT_DIR)/../config.mk

CFLAGS  := -c -fPIC -std=c++0x -m64

ifeq ($(DEBUG),)
CFLAGS  := $(CFLAGS) -O2 -DNDEBUG
else
CFLAGS  := $(CFLAGS) -g
endif

COMPILE := $(CFLAGS) \
	   -I$(ROOT_DIR)/include/kflow \
	   -I$(BOOST_DIR)/include \
	   -I$(GFLAGS_DIR)/include \
	   -I$(GLOG_DIR)/include

DEPS    := src/Common.o \
	   src/Pipeline.o \
           src/Stage.o 

DEST	:= $(ROOT_DIR)/lib/libkflow.a

all: $(DEST)

debug:
	$(MAKE) DEBUG=1 all

$(DEST): $(DEPS)
	$(AR) rvs $@ $(DEPS)

%.o: %.cpp $(ROOT_DIR)/include/kflow/*.h
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf $(DEPS)
	rm -rf $(DEST)

.PHONY: all debug clean

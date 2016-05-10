ROOT_DIR:=..
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

LINK	:= -L$(BOOST_DIR)/lib \
           -lboost_system -lboost_thread -lboost_iostreams -lboost_filesystem -lboost_regex \
	   -L$(GFLAGS_DIR)/lib -lgflags \
           -L$(GLOG_DIR)/lib -lglog \
           -lpthread -lm -ldl

DEPS    := Common.o \
	   Pipeline.o \
           Stage.o 

DEST	:= $(ROOT_DIR)/lib/libkflow.so

all: $(DEST)

debug:
	$(MAKE) DEBUG=1 all

$(DEST): $(DEPS)
	$(PP) -shared -o $@ $(DEPS) $(LINK) 

%.o: %.cpp $(ROOT_DIR)/include/kflow/*.h
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf *.o
	rm -rf $(DEST)

.PHONY: all debug clean

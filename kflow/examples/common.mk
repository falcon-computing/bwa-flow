ROOTDIR := ../../..
include $(ROOTDIR)/config.mk

# check all variables
ifeq ($(BOOST_DIR),)
$(error BOOST_DIR is not set properly in Makefile.config)
endif
ifeq ($(GLOG_DIR),)
$(error GLOG_DIR is not set properly in Makefile.config)
endif

# binaries
PP=g++
CC=gcc
CC=ar
RM=rm
ECHO=echo
MAKE=make

CFLAGS   := -c -fPIC -std=c++0x $(CFLAGS)

ifeq ($(DEBUG),)
CFLAGS   := $(CFLAGS) -O2 -DNDEBUG
else
CFLAGS   := $(CFLAGS) -g
endif

INCLUDES := -I$(ROOTDIR)/kflow/include \
	    -I$(BOOST_DIR)/include \
	    -I$(GLOG_DIR)/include \
	    -I$(GFLAGS_DIR)/include \
	    $(INCLUDES)

COMPILE := $(CFLAGS) $(INCLUDES)

LINK    := -L$(ROOTDIR)/kflow/lib -lkflow \
	   -L$(BOOST_DIR)/lib \
	       -lboost_system \
	       -lboost_thread \
	       -lboost_iostreams \
	       -lboost_filesystem \
	       -lboost_regex \
	   -L$(GLOG_DIR)/lib -lglog \
	   -L$(GFLAGS_DIR)/lib -lgflags \
	   -lpthread -lm -ldl \
	   $(LINK)

all: $(DST)

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	$(RM) -rf *.o
	$(RM) -rf $(DST)

.PHONY: all clean

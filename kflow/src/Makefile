ROOT_DIR=..
include $(ROOT_DIR)/../config.mk

CFLAGS  := -c -fPIC -std=c++0x -m64 -g 
COMPILE := $(CFLAGS) \
				 	 -I$(ROOT_DIR)/include \
				 	 -I$(BOOST_DIR)/include \
				 	 -I$(GLOG_DIR)/include
LINK		:= -L$(BOOST_DIR)/lib \
					 -lboost_system -lboost_thread -lboost_iostreams -lboost_filesystem -lboost_regex \
					 -L$(GLOG_DIR)/lib -lglog \
					 -lpthread -lm -ldl

DEPS		:= Pipeline.o \
					 Stage.o 

DEST		:= $(ROOT_DIR)/lib/libkflow.so

all: $(DEST)

$(DEST): $(DEPS)
	$(PP) -shared -o $@ $(DEPS) $(LINK) 

%.o: %.cpp $(ROOT_DIR)/include/*.h
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf *.o
	rm -rf $(DEST)

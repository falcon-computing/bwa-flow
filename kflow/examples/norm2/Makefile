ROOT_DIR=../..
include $(ROOT_DIR)/../config.mk

CFLAGS  := -c -fPIC -std=c++0x -m64 -g 
COMPILE := $(CFLAGS) \
				 	 -I$(ROOT_DIR)/include \
				 	 -I$(BOOST_DIR)/include \
				 	 -I$(GLOG_DIR)/include
LINK		:= -L$(ROOT_DIR)/lib -lkflow \
				   -L$(BOOST_DIR)/lib \
					 -lboost_system -lboost_thread -lboost_iostreams -lboost_filesystem -lboost_regex \
					 -L$(GLOG_DIR)/lib -lglog \
					 -lpthread -lm -ldl

all: compute gen

compute: compute.o 
	$(PP) -o $@ $< $(LINK) 

gen: gen.o 
	$(PP) -o $@ $< $(LINK) 

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf *.o
	rm -rf app 

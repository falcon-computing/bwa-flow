ROOT_DIR=../..
include $(ROOT_DIR)/make.config

CFLAGS  := -shared -fPIC -std=c++0x -m64 -g 
COMPILE := $(CFLAGS) \
				 	 -I$(ROOT_DIR)/include \
				 	 -I$(BOOST_DIR)/include \
				 	 -I$(GLOG_DIR)/include
LINK		:= -L$(BOOST_DIR)/lib \
					 -lboost_system -lboost_thread -lboost_iostreams -lboost_filesystem -lboost_regex \
					 -L$(GLOG_DIR)/lib -lglog \
					 -lpthread -lm -ldl

all: app

app: main.o 
	$(PP) -o $@ main.o $(LINK) 

%.o: %.cpp Makefile
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf *.o
	rm -rf app 

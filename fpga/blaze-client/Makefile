include /curr/diwu/prog/blaze/Makefile.config

MANAGER_DIR=/curr/diwu/prog/blaze/manager

PP	:= g++

CFLAGS	:= -s -fPIC -std=c++0x -g
COMPILE	:= -c $(CFLAGS) \
	   -I$(MANAGER_DIR)/include \
	   -I$(BOOST_DIR)/include \
	   -I$(PROTOBUF_DIR)/include \
	   -I$(GLOG_DIR)/include \
	   -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux 
	
LINK	:= -L$(MANAGER_DIR)/lib -lblaze \
	   -L$(BOOST_DIR)/lib \
	   	-lboost_system \
		-lboost_thread \
		-lboost_iostreams \
		-lboost_filesystem \
		-lboost_regex \
	   -L$(PROTOBUF_DIR)/lib -lprotobuf \
	   -L$(GLOG_DIR)/lib -lglog \
	   -lpthread -lm -ldl

all: app

app: main.o 
	$(PP) -o $@ main.o $(LINK) 

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf *.o
	rm -rf app 

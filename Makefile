BWA_DIR := ./bwa
SRC_DIR := ./src
include /curr/diwu/prog/blaze/Makefile.config

MANAGER_DIR=/curr/diwu/prog/blaze/manager


CC	:= gcc
PP	:= g++

CFLAGS 	:= -g -std=c++0x -fPIC -O3
OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
           $(SRC_DIR)/chain2alnhw.o\
           $(SRC_DIR)/smithwaterman.o\
	   $(SRC_DIR)/main.o \
	   $(SRC_DIR)/SWRead.o \
	   $(SRC_DIR)/util.o

PROG	:= ./bin/bwa
INCLUDES:= -I. -I./bwa -I$(MANAGER_DIR)/include \
	   -I$(BOOST_DIR)/include \
           -I/curr/software/Xilinx/Vivado_HLS/2015.4/include\
	   -I$(PROTOBUF_DIR)/include \
	   -I$(GLOG_DIR)/include \
	   -I$(JAVA_HOME)/include -I$(JAVA_HOME)/include/linux 
	
LIBS	:= -L$(BWA_DIR) -lbwa \
	   -L$(MANAGER_DIR)/lib -lblaze \
	   -L$(BOOST_DIR)/lib \
	   	-lboost_system \
		-lboost_thread \
		-lboost_iostreams \
		-lboost_filesystem \
		-lboost_regex \
	   -L$(PROTOBUF_DIR)/lib -lprotobuf \
	   -L$(GLOG_DIR)/lib -lglog \
	   -lpthread -lm -ldl -lz -lrt



all:$(PROG)

./bin/bwa: $(BWA_DIR)/libbwa.a $(OBJS) 
	$(PP) $(OBJS) -o $@ $(LIBS)

$(SRC_DIR)/%.o:	$(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@
$(SRC_DIR)/%.o:	$(SRC_DIR)/%.cpp
	$(PP) -c $(CFLAGS) $(INCLUDES) $< -o $@

./bwa/libbwa.a:
	make -C $(BWA_DIR)

clean:
	rm -f gmon.out 
	rm -f $(OBJS)
	rm -f $(PROG)  

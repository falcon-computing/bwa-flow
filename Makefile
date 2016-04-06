BWA_DIR := ./bwa
SRC_DIR := ./src

CC	:= gcc
PP	:= g++

CFLAGS 	:= -g -std=c++0x -fPIC -O2
OBJS	:= $(SRC_DIR)/wrappered_mem.o \
	   $(SRC_DIR)/preprocess.o \
	   $(SRC_DIR)/main.o \
	   $(SRC_DIR)/util.o

PROG	:= ./bin/bwa
INCLUDES:= -I.
LIBS	:= -L$(BWA_DIR) -lbwa \
	   -lm -lz -lpthread -lrt

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

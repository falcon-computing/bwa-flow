include config.mk

CFLAGS		:= -Wall -Wno-unused-function -O2
WRAP_MALLOC	:= -DUSE_MALLOC_WRAPPERS
DFLAGS		:= -DHAVE_PTHREAD $(WRAP_MALLOC)

INCLUDES	:= -Ibwa
LIBS		:= -lm -lz -lpthread

COMPILE		:= -c $(CFLAGS) $(INCLUDES)
DST			:= ./bin/bwa
OBJS		:= ./src/bntseq_newflow.o \
			   ./src/bwa_newflow.o \
			   ./src/bwt_newflow.o \
			   ./src/kopen.o \
			   ./src/kstring.o \
			   ./src/ksw_newflow.o \
			   ./src/main.o \
			   ./src/mem_extra_newflow.o \
			   ./src/mem_newflow.o \
			   ./src/mem_pair_newflow.o \
			   ./src/preprocess.o \
			   ./src/utils.o

all: $(DST)

./bin/bwa: $(OBJS)
	$(PP) $(CFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)

%.o: %.cpp
	$(PP) $(COMPILE) $< -o $@

clean:
	rm -rf $(OBJS)
	rm -rf $(DST)

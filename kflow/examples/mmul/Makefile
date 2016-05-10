DST=app \
    base
			
MKLROOT	:= /curr/diwu/intel/mkl
INCLUDES:= -I$(MKLROOT)/include
LINK	:= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
#CFLAGS	:= -p

include ../common.mk

app: main.o 
	$(PP) -p -o $@ $< $(LINK) 

base: base.o 
	$(PP) -p -o $@ $< $(LINK) 

#__INSERT_LICENSE__
.PHONY: all

default: run_nettest

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 


#w 
femref.bin: femref.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)


#s

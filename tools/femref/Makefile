#__INSERT_LICENSE__
.PHONY: all

default: run_nettest

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

SRCS := femref.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w 
femref.bin: femref.o
	${CXX_CLINKER} -o $@ $<  $(LDFLAGS_NO_PROG_LIB)

hasher.o: hasher.cpp	
	g++ -g -c -o $@ $<

myhash.o: myhash.cpp	
	g++ -g -c -o $@ $<

myhash.bin: myhash.o hasher.o
	g++ -o $@ $^

#s

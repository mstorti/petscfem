#__INSERT_LICENSE__
.PHONY: all

default: run_nettest

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

SRCS := femref.cpp hasher.cpp myhash.cpp myhash2.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 
include makefilec.d

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

#w 
femref.bin: femref.o # hasher.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

hasher.o: hasher.cpp	
	g++ -g -c -o $@ $<

myhash.o: myhash.cpp	
	g++ -g -c -o $@ $<

myhash.bin: myhash.o hasher.o
	g++ -o $@ $^

myhash2.bin: myhash2.o hasher.o
	g++ -o $@ $^

myhash3.bin: myhash3.o hasher.o
	g++ -o $@ $^

md5.bin:
	gcc -o md5 md5.c md5r.c

trymd5.bin: trymd5.c md5r.c
	gcc -g -o $@ $^

cdepend:
	gcc -MM *.c > makefilec.d

#s

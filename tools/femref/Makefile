#__INSERT_LICENSE__
.PHONY: all

N := 5
OPT := -O2 -funroll-loops
# OPT := -g

PROG := myhash$(N).bin
default: $(PROG)

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

SRCS := femref.cpp femref2.cpp hasher.cpp myhash.cpp myhash2.cpp	\
	myhash3.cpp myhash4.cpp myhash5.cpp gtemplates.cpp 		\
	visitor.cpp combiner.cpp trytree.cpp

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 
include makefilec.d

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

trytree.bin: trytree.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

femref.bin: femref.o femref2.o hasher.o md5r.o gtemplates.o \
		visitor.o combiner.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

hasher.o: hasher.cpp
	g++ $(OPT) -c -o $@ $<

myhash.o: myhash.cpp	
	g++ $(OPT) -c -o $@ $<

myhash.bin: myhash.o hasher.o
	g++ -o $@ $^

myhash2.bin: myhash2.o hasher.o
	g++ -o $@ $^

md5r.o: md5r.c
	gcc $(OPT) -c -o $@ $<

myhash$(N).bin: myhash$(N).o hasher.o md5r.o
	g++ $(OPT) -o $@ $^

md5.bin:
	gcc -o md5 md5.c md5r.c

trymd5.bin: trymd5.c md5r.c
	gcc $(OPT) -o $@ $^

cdepend:
	gcc -MM *.c > makefilec.d

#s

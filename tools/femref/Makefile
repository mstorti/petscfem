#__INSERT_LICENSE__
.PHONY: all

N := 5
OPT := -O2 -funroll-loops
# OPT := -g

PROG := myhash$(N).bin
default: getsurf.bin

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

CPP_FILES := femref2 hasher md5r gtemplates visitor combiner

SRCS := $(CPP_FILES=.cpp) getsurf.cpp femref.cpp  myhash.cpp myhash2.cpp \
	myhash3.cpp myhash4.cpp myhash5.cpp trytree.cpp 

BOPT := O_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 
include makefilec.d

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

tryme2.bin: tryme2.o
	g++ -o $@ $^

trytree.bin: trytree.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

femref.bin: $(CPP_FILES:=.o)
	echo ${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

getsurf.bin: getsurf.o $(CPP_FILES:=.o)
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

getsurf:
	$(GETSURF) -b 1 -x $(DATA)/cube.nod.tmp	\
	-c $(DATA)/proy.con.tmp -u ./cube.state.tmp
	-s proy.surf-con.tmp -g proy.grad-u.tmp

#s

#__INSERT_LICENSE__
.PHONY: all

N := 5
OPT := -O2 -funroll-loops
# OPT := -g

PROG := myhash$(N).bin
default: libfemref.so
# default: getsurf.bin

scripts = odoc.pl makewhat.pl

#v $hdrmxlen += 3;
#p [in tools]
#w Builds all doc in this directory
all: $(scripts:pl=1) $(scripts:pl=txt) 

dvectori.cpp: dvector-epl.cpp dvector.cpp
	make-depl $< $@ -d type=int -p '^\s*//PP>'

dvectord.x: dvectord.cpp

dvectord.cpp: dvector-epl.cpp dvector.cpp
	make-depl $< $@ -d type=double -d sfx=dbl -p '^\s*//PP>'

dvectori.o: dvectori.cpp dvectori.x
dvectord.o: dvectord.cpp dvectord.x

CPP_FILES := femref2 hasher md5r gtemplates visitor combiner \
		dvectori dvectord dvector2 gsguile

SRCS := $(CPP_FILES=.cpp) getsurf.cpp femref.cpp  myhash.cpp myhash2.cpp \
	myhash3.cpp myhash4.cpp myhash5.cpp trytree.cpp 

BOPT := g_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 
CPPFLAGS += -Wno-deprecated

ifeq (makefilec.d,$(wildcard makefilec.d)) 
include makefilec.d
endif

ifeq (makefile.d,$(wildcard makefile.d)) 
include makefile.d
endif

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

%.x: %.cpp
	guile-snarf -o $@ $(CPPFLAGS) $<

gsguile.o: gsguile.x

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

bjhash.o: bjhash.c
	gcc -c -o $@ $<
myhash$(N).bin: myhash$(N).o hasher.o md5r.o bjhash.o
	g++ $(OPT) -o $@ $^

md5.bin:
	gcc -o md5 md5.c md5r.c

trymd5.bin: trymd5.c md5r.c
	gcc $(OPT) -o $@ $^

cdepend:
	gcc -MM *.c > makefilec.d

GETSURF := ./getsurf.bin
cube:
	make-depl cube.dx.epl cube2.dx

getsurf:
ifeq (1,0)
	$(GETSURF) -b 1 -x ./cube.nod.tmp	\
	-c ./cube.con.tmp -u ./cube.state.tmp	\
	-s cube.surf-con.tmp -g cube.grad-u.tmp
else
	guile -s femref.scm
endif

ifeq (1,1)
# LD_LIBRARY_PATH += ../../src/libpetscfem_g.a
libfemref.so: getsurf.o $(CPP_FILES:=.o)
	${CXX_CLINKER} -g -shared				\
		-o $@ getsurf.o $(CPP_FILES:=.o)		\
		`guile-config link` $(LDFLAGS_NO_PROG_LIB)
else
libfemref.so: getsurf2.o
	${CXX_CLINKER} -g -shared			\
		-o $@ getsurf2.o `guile-config link`	\
		$(LDFLAGS_NO_PROG_LIB)
endif

LDFLAGS := 						\
	-lstdc++ -Wl,-E,-rpath,$(LD_LIBRARY_PATH)	\
	-L$(PETSCFEM_DIR)/src				\
	-lpetscfem$(osfx)				\
	-L$(NEWMAT)  -lnewmat				\
	-lglib						\
	-L$(LIBRETTO_LIB) -libretto			\
	-L$(PETSC_DIR)/lib/lib$(BOPT)/linux/		\
	$(PETSC_SLES_LIB)				\
	$(SUPERLU)/libsuperlu$(osfx).a			\
	$(MESCHACH)/libmes.a				\
	-L$(METIS) -lmetis				\
	$(SSL)/simpleskts.a				\
	-lpthread					\
	-lblas -llapack 

mpiguile: mpiguile.o
	gcc -g -o $@ $< `guile-config link` 	\
	$(LDFLAGS)

mpitest: mpitest.o
	gcc -g -o $@ $< `guile-config link` 	\
	$(LDFLAGS)

yuv_compress:
	@echo nothing

#s

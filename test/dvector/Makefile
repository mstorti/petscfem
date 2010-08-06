#__INSERT_LICENSE__
.PHONY: all

default: tests

tests: test_dvclonepar

SRCS := dvectest.cpp

# BOPT := g_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

ifeq (makefile.d,$(wildcard makefile.d)) 
include makefile.d
endif

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

mylib:
	$(MAKE) -C $(PETSCFEM_DIR)/src 

prog: dvectest.bin

reprog: 
	$(MAKE) -W dvectest.cpp dvectest.bin

dvectest.bin: dvectest.o
	$(MAKE) mylib
	${CXX_CLINKER} -o $@ $^ $(LDFLAGS_NO_PROG_LIB)

local_clean::
	@rm -f dvectest.bin dvectest.o dvectest.out

test_dvclonepar:
	@rm -f dvectest.out
	$(MAKE) reprog
	$(MPIEXEC) $(MCH) -n 4 dvectest.bin > dvectest.out

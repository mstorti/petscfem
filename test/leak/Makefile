#__INSERT_LICENSE__
# $Id: Makefile,v 1.3 2002/10/07 00:25:32 mstorti Exp $ 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ifeq (1,1)
CPPFLAGS += -DUSE_OLD_PETSC_VERSION
PETSC_DIR = $(HOME)/PETSC/petsc-2.0.24
include ${PETSC_DIR}/bmake/$(PETSC_ARCH)/base
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
else
PETSC_DIR = $(HOME)/PETSC/petsc-2.1.3
include ${PETSC_DIR}/bmake/common/base
endif

BOPT=g_c++
LD_LIBRARY_PATH=$(PETSC_DIR)/lib/lib$(BOPT)/$(PETSC_ARCH)
LDFLAGS = -Wl,-E,-rpath,$(LD_LIBRARY_PATH) ${PETSC_SLES_LIB} -lc -lm -lg2c

leak.bin: leak.o
	-rm -f $@
	${CXX_CLINKER} -o $@ leak.o $(LDFLAGS)

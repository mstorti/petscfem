 # $Id: Makefile,v 1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

MYOBJS = genload.o lapla.o laplace.o
SRCS := $(MYOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

compile: $(PROG).bin

lclean: general_clean

# File containing dependencies computed by makedepend
include makefile.d

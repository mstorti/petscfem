 # $Id: Makefile,v 1.2 2001/01/02 16:50:30 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

MYOBJS = genload.o lapla.o laplace.o
SRCS := $(MYOBJS:.o=.cpp)

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

include $(PETSCFEM_DIR)/Makefile.base

compile: $(PROG).bin

lclean: general_clean

# File containing dependencies computed by makedepend
#

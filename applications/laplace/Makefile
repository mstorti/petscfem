#__INSERT_LICENSE__
 # $Id: Makefile,v 1.4 2001/04/01 03:10:44 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

MYOBJS = genload.o lapla.o laplace.o
SRCS := $(MYOBJS:.o=.cpp)

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

include $(PETSCFEM_DIR)/Makefile.base

compile: $(PROG).bin

#lclean: general_clean

# File containing dependencies computed by makedepend
#

#__INSERT_LICENSE__
 # $Id: Makefile,v 1.6 2002/09/17 17:15:13 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

LIBOBJS = genload.o lapla.o laplace.o compdate.o
MYOBJS = 
SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

#__INSERT_LICENSE__
 # $Id: Makefile,v 1.5 2002/07/25 22:35:31 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

LIBOBJS = genload.o lapla.o laplace.o 
MYOBJS = 
SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

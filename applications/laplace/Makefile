#__INSERT_LICENSE__
 # $Id: Makefile,v 1.6 2002/09/17 17:15:13 mstorti Exp $ 

.PHONY: all run lclean save compile

PROG = laplace
PETSCFEM_DIR   = ../..

PROGOBJ = main.o
LIBOBJS = laplace.o compdate.o genload.o lapla.o
MYOBJS = 
SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

TAGFLAGS = $(TAGINCL)="../../src/TAGS"

SRCS := $(MYOBJS:.o=.cpp) $(LIBOBJS:.o=.cpp)

include $(PETSCFEM_DIR)/Makefile.base

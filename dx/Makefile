#-*- mode: Makefile -*-
SHELL = /bin/sh
BASE = /usr/local/dx

default: epimport

# need arch set, e.g. by
# setenv DXARCH `dx -whicharch`
include $(BASE)/lib_$(DXARCH)/arch.mak
DX_RTL_LDFLAGS := -shared
OBJEXT := o

PETSCFEM_DIR := ..
include $(PETSCFEM_DIR)/Makefile.base

FILES_epimport = userepimport.$(OBJEXT) epimport.$(OBJEXT)

BIN = $(BASE)/bin

CFLAGS = -I./ -I$(BASE)/include $(DX_CFLAGS) -DUSE_SSL -I$(SSL) -I$(PETSCFEM_DIR)

LDFLAGS = -L$(BASE)/lib_$(DXARCH)

LIBS = -lDX $(DX_GL_LINK_LIBS) $(DXEXECLINKLIBS) 

OLIBS = -lDXlite -lm

BIN = $(BASE)/bin

# create the necessary executable
epimport: $(FILES_epimport) 
	$(SHARED_LINK) $(DXABI) $(LDFLAGS) -o epimport userepimport.$(OBJEXT)	\
		epimport.$(OBJEXT) $(DX_RTL_LDFLAGS) $(SYSLIBS)			\
		$(SSL)/simpleskts.a $(LIBPETSCFEM) -lstdc++

epimport.o: epimport.cpp
	g++ -c $(DXABI) $(DX_RTL_CFLAGS) $(CFLAGS) epimport.cpp

# a command to run the user module
run: epimport 
	dx -edit -mdf epimport.mdf &

# make the user files
userepimport.c: epimport.mdf
	$(BIN)/mdf2c -m epimport.mdf > userepimport.c

userepimport.o: userepimport.c
	cc -c $(DXABI) $(DX_RTL_CFLAGS) $(CFLAGS) userepimport.c 

local_clean::
	-rm -f epimport

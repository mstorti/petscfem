#-*- mode: Makefile -*-
SHELL = /bin/sh
BASE = /usr/local/dx

# need arch set, e.g. by
# setenv DXARCH `dx -whicharch`
include $(BASE)/lib_$(DXARCH)/arch.mak
DX_RTL_LDFLAGS := -shared
OBJEXT := o
SSL := $(HOME)/SOFT/SSL

FILES_epimport = userepimport.$(OBJEXT) epimport.$(OBJEXT)

BIN = $(BASE)/bin

CFLAGS = -I./ -I$(BASE)/include $(DX_CFLAGS) -I$(SSL)

LDFLAGS = -L$(BASE)/lib_$(DXARCH)

LIBS = -lDX $(DX_GL_LINK_LIBS) $(DXEXECLINKLIBS) 

OLIBS = -lDXlite -lm

BIN = $(BASE)/bin

# create the necessary executable
epimport: $(FILES_epimport) 
	$(SHARED_LINK) $(DXABI) $(LDFLAGS) -o epimport userepimport.$(OBJEXT) \
		epimport.$(OBJEXT) $(DX_RTL_LDFLAGS) $(SYSLIBS)  $(SSL)/simpleskts.a

.c.o: ; cc -c $(DXABI) $(DX_RTL_CFLAGS) $(CFLAGS) $*.c 

.C.o: ; cc -c $(DXABI) $(DX_RTL_CFLAGS) $(CFLAGS) $*.C 

# a command to run the user module
run: epimport 
	dx -edit -mdf epimport.mdf &

# make the user files
userepimport.c: epimport.mdf
	$(BIN)/mdf2c -m epimport.mdf > userepimport.c

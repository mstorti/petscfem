# $Id: Makefile,v 1.32 2003/03/23 20:52:43 mstorti Exp $ 

.PHONY: force all run tests default strip gasflow nozzle  vtubei

default: vtube

tests: 

#w Gas flow example (basic)
gasflow:
	$(MAKE) CASE_NAME=gasflow_test_2d run

#w Nozzle example
nozzle:
	$(MAKE) CASE_NAME=$@ run

#w Vortex tube example
vtube: gasflow.efn
	$(MAKE) CASE_NAME=$@ run

#w Vortex tube example (incompressible)
vtubei:
	$(MAKE) CASE_NAME=vtube incomp=1 PROG=$(NS_O) run

#w Strip example
strip:
	$(MAKE) CASE_NAME=$@ run

SRCS = 
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
include $(PETSCFEM_DIR)/test/Makefile.base.test

ifeq ($(CASE_NAME),) 
CASE_NAME := strip
endif
LOC_PROG := advdif
ifeq ($(PROG),)
PROG := $(ADVDIF)
endif

nolocal := yes
include ../les/stabi3d/Makefile.parproc
# NP := 2
run: $(CASE_NAME).depl
	-$(MPIRUN) $(NL) -np $(NP) -machinefile \
		machi.dat $(PROG) -case $< 2>&1

$(CASE_NAME).depl: force

local_clean::
	-rm -rf mat.output visc2inv inv2visc

#w generate DX files
dx:
	$(EPERL) -P viewer.dx > viewer_tmp.dx

#w Make extension function (hook and dynamic temporal function)
gasflow.efn: gasflow.cpp

# move states
move:
	find .  -maxdepth 1 -name 'vtube.state_*.tmp' | xargs -n 1 echo gzip
	stmove.pl -f ./ -t RUN -p 'vtube.state_*.tmp.gz'

epi:
	$(MAKE) -f epimport.make epimport

socket.o: socket.cpp
socket.bin: socket.o
	$(CXX) -o $@ $<

SSL := $(HOME)/SOFT/SSL
socket2.o: CPPFLAGS = -I$(SSL)
socket2.o: socket2.cpp
socket2.bin: LDFLAGS = $(SSL)/simpleskts.a -lpthread
socket2.bin: socket2.o
	$(CXX) -o $@ $< $(LDFLAGS)


#s

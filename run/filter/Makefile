PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base


LDLIBS = -lm

CPPFLAGS = -g

LDFLAGS = -g

run: sttfilter
	sttfilter

sttfilter: sttfilter.o
	g++ $(LDFLAGS) -o $@  $^ $(LDLIBS)

sttfilter.o: sttfilter.h sttfilter.cpp

local_clean::
	-rm -f filter.out sttfilter

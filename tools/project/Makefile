#__INSERT_LICENSE__
.PHONY: all

default: project3.bin

SRCS := project.cpp project2.cpp project3.cpp project4.cpp \
	project-wrp.cpp project5.cpp

# BOPT := O_c++
PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

ifeq (makefile.d,$(wildcard makefile.d)) 
include makefile.d
endif

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

project.bin: project.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

project2.bin: project2.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

project3.bin: project3.o project.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

#w Make "project4.bin" which projects results from a FEM mesh
# to another
project4.bin: project4.o project.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

project5.bin: project5.o project.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

yuv_compress:
	@echo nothing

project-wrp.x: project-wrp.cpp

libproject.so: project.o project-wrp.o
	${CXX_CLINKER} -g -shared \
		-o $@ $^ `guile-config link` $(LDFLAGS_NO_PROG_LIB)

fluent.dx: fluent-in.dx
	make-depl $< $@

#s

what:
	@makewhat.pl -s Makefile

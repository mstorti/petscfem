#__INSERT_LICENSE__
.PHONY: all

default: project3.o

SRCS := project.cpp project2.cpp project3.cpp project4.cpp

BOPT := O_c++
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

project3.bin: project3.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

project4.bin: project4.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

yuv_compress:
	@echo nothing

#s

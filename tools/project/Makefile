#__INSERT_LICENSE__
.PHONY: all

OPT := -O2 -funroll-loops

default: project.bin

BOPT := g_c++
include $(PETSCFEM_DIR)/Makefile.base
OMKFLS = 

ifeq (makefile.d,$(wildcard makefile.d)) 
include makefile.d
endif

CLEAN_DIRS = 
DEPEND_DIRS = $(CLEAN_DIRS)

project.bin: project.o
	${CXX_CLINKER} -o $@ $^  $(LDFLAGS_NO_PROG_LIB)

yuv_compress:
	@echo nothing

#s

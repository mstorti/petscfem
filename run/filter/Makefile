#__INSERT_LICENSE__
run: sttfilter
	sttfilter

PETSCFEM_DIR = ../..
include $(PETSCFEM_DIR)/Makefile.base


LDLIBS = -lm

CPPFLAGS = -g

LDFLAGS = -g

SRCDIR = $(PETSCFEM_DIR)/src

STTFLT = libpetscfem.a(sttfilter.o)

$(SRCDIR)/$(STTFLT):
	$(MAKE) -C $(SRCDIR) $(STTFLT)

sttfilter: sttfilter.o $(SRCDIR)/libpetscfem.a
	g++ $(LDFLAGS) -o $@  $^ $(LDLIBS)

sttfilter.o: $(SRCDIR)/sttfilter.h

local_clean::
	-rm -f filter.out sttfilter

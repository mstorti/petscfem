 # $Id: Makefile,v 1.1.1.1 2000/12/28 12:54:43 mstorti Exp $ 

.PHONY: run compile static

BOPT = g

PROG = petscfem
PETSCFEM_DIR   = ..
MYOBJS = petscfem.o petscfem_wrap.o

#include $(PETSCFEM_DIR)/Makefile.base

include $(PETSCFEM_DIR)/../Makefile.defs

CPPFLAGS         = $(RH60) ${BS_INCLUDE} -I$(GLIB) -I$(NEWMAT) \
                     -I$(METIS)/Lib \
                    -I$(PETSCFEM_DIR) -I$(MESCHACH) \
                    -I$(LIBRETTO_I)
FPPFLAGS         =

include ${PETSC_DIR}/bmake/${PETSC_ARCH}/base

CFLAGS = -Dbool=char -DHAS_BOOL -I/usr/lib/perl5/5.00503/i386-linux/CORE 

# El comando que anduvo con el ejemplo es
# gcc -Dbool=char -DHAS_BOOL -I/usr/lib/perl5/5.00503/i386-linux/CORE \
#               -c example.c example_wrap.c 

PETSC_STATIC_LIB =							\
-L./petsc -lpetscsles	\
 -lpetscdm								\
 -lpetscmat								\
 -lpetscvec								\
 -lpetsc /usr/local/lib/liblapack.so.2.0.1				\
/usr/local/lib/libblas.so.2.0.1 -L/usr/local/mpich-1.2.0/lib -lpmpich	\
-lmpich -L/usr/X11R6/lib -lX11 -lc -lg2c -lm -lm

$(PROG)_wrap.c: $(PROG).i
	swig -static -perl5 $(PROG).i

$(PROG).so: $(PROG).o $(PROG)_wrap.o 
	ld -G $(PROG).o $(PROG)_wrap.o -o $(PROG).so -lc \
		-L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/ \
		-L/home/mstorti/PETSC/petsc-2.0.24/lib/libg/linux/ \

static: $(PROG).o $(PROG)_wrap.o
	 ld $(PROG).o $(PROG)_wrap.o -o myperl -lc \
		-L/usr/lib/perl5/5.00503/i386-linux/CORE/ \
              -lperl -lm -o myperl \
		-L/usr/lib/gcc-lib/i386-redhat-linux/egcs-2.91.66/ \
		-lperl -lm ${PETSC_STATIC_LIB} 

compile: $(PROG).bin

$(PROG).bin: $(MYOBJS) 
	-rm $(PROG).bin
	${C_CLINKER} -o $(PROG).bin $(MYOBJS) -lc ${PETSC_SLES_LIB} 

run: $(PROG).bin
	-${MPIRUN} -np 1 $(ADV)

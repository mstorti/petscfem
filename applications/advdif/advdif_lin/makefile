# $Id: makefile,v 1.1.2.1 2004/03/01 21:55:01 mstorti Exp $ 

# -----------------------------------------------------------------------
# Specify  user-defined compiler flags
# -----------------------------------------------------------------------
CFLAGS    =
# -----------------------------------------------------------------------
# Define application source and examples
# -----------------------------------------------------------------------
#

ALL: main

LOCDIR         = /u/rodrigop/petscfem-hydro/applications/advdif/advdif_lin/
EXAMPLESC      = main.c appview.c appload.c appalgebra.c appelement.c
EXAMPLESO      = main.o appview.o appload.o appalgebra.o appelement.o
EXAMPLECH      = appctx.h
DOCS           = finput.c
DIRS           = scalable

advdif_lin_nn.bin: chkopts ${EXAMPLESO}
	-$(RM) advdif_lin_nn.bin
	-$(CLINKER) -o advdif_lin_nn.bin ${EXAMPLESO} $(PETSC_SNES_LIB)
	-$(RM) ${EXAMPLESO}

finput: chkopts finput.o
	-$(RM) finput
	-$(CLINKER) -o finput finput.o $(PETSC_SNES_LIB)
	-$(RM) finput.o


include $(PETSC_DIR)/bmake/common/base

TESTEXAMPLES_1 = main.PETSc main.rm

include $(PETSC_DIR)/bmake/common/test



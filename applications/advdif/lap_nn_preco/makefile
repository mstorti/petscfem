# $Id: makefile,v 1.1.2.1 2004/03/02 17:02:24 mstorti Exp $ 

# -----------------------------------------------------------------------
# Specify  user-defined compiler flags
# -----------------------------------------------------------------------
CFLAGS    =
# -----------------------------------------------------------------------
# Define application source and examples
# -----------------------------------------------------------------------
#

ALL: main

LOCDIR         = /u/rodrigop/petsfem-hydro/applications/advdif_nnpreco/
EXAMPLESC      = main.c appview.c appload.c appalgebra.c appelement.c
EXAMPLESO      = main.o appview.o appload.o appalgebra.o appelement.o
EXAMPLECH      = appctx.h
DOCS           = finput.c
DIRS           = scalable

lap_nn.bin: chkopts ${EXAMPLESO}
	-$(RM) lap_nn.bin
	-$(CLINKER) -o lap_nn.bin ${EXAMPLESO} $(PETSC_SNES_LIB)
	-$(RM) ${EXAMPLESO}

finput: chkopts finput.o
	-$(RM) finput
	-$(CLINKER) -o finput finput.o $(PETSC_SNES_LIB)
	-$(RM) finput.o


include $(PETSC_DIR)/bmake/common/base

TESTEXAMPLES_1 = main.PETSc main.rm

include $(PETSC_DIR)/bmake/common/test



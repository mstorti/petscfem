//__INSERT_LICENSE__
//$Id: compdate.cpp,v 1.4.100.1 2007/02/19 20:23:56 mstorti Exp $
#include <cstdio>
#include <petsc.h>
#include "./pfversion.h"

extern MPI_Comm PETSCFEM_COMM_WORLD;

#ifndef PETSCFEM_HOSTNAME
#define PETSCFEM_HOSTNAME "<unknown host>"
#endif

void petscfem_print_link_date() {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "-----------------------------------------------------------------\n"
	      "Navier-Stokes module built on " __DATE__ ", " __TIME__ ", with "
#ifdef __GNUC__
	      "GCC %d.%d.%d\n"
#else
	      "<unknown compiler>\n"
#endif
	      "on host: " PETSCFEM_HOSTNAME "\n",
#ifdef __GNUC__
	      __GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__
#endif
	      );
  PetscPrintf(PETSCFEM_COMM_WORLD,
              "PETSc-FEM Version %s (Git)\n Full Commit: %s\nDate: %s\n"
	      "-----------------------------------------------------------------\n",
              DOCVERSION,DOCCOMMIT,DOCDATE);
}

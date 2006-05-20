//__INSERT_LICENSE__
//$Id: compdate.cpp,v 1.4.70.2 2006/05/20 21:44:40 dalcinl Exp $
#include <cstdio>
#include <petsc.h>

#ifndef PETSCFEM_HOSTNAME
#define PETSCFEM_HOSTNAME "<unknown host>"
#endif

extern MPI_Comm PETSCFEM_COMM_WORLD;

void petscfem_print_link_date() {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "-----------------------------------------------------------------\n"
	      "Navier-Stokes module built on " __DATE__ ", " __TIME__ ", with "
#ifdef __GNUC__
	      "GCC %d.%d.%d\n"
#else
	      "<unknown compiler>\n"
#endif
	      "on host: " PETSCFEM_HOSTNAME "\n"
	      "-----------------------------------------------------------------\n",
#ifdef __GNUC__
	      __GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__
#endif
	      );
}

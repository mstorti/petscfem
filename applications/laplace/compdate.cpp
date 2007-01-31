//__INSERT_LICENSE__
//$Id: compdate.cpp,v 1.2.92.1 2007/01/31 02:02:56 dalcinl Exp $
#include <cstdio>
#include <petsc.h>

extern MPI_Comm PETSCFEM_COMM_WORLD;

#ifndef PETSCFEM_HOSTNAME
#define PETSCFEM_HOSTNAME "<unknown host>"
#endif

void petscfem_print_link_date() {
  PetscPrintf(PETSCFEM_COMM_WORLD,
	      "-----------------------------------------------------------------\n"
	      " module built on " __DATE__ ", " __TIME__ ", with "
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

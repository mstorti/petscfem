//__INSERT_LICENSE__
//$Id: compdate.cpp,v 1.1 2002/09/17 15:25:02 mstorti Exp $
#include <cstdio>
#include <petsc.h>

#ifndef PETSCFEM_HOSTNAME
#define PETSCFEM_HOSTNAME "<unknown host>"
#endif

void print_petscfem_link_date() {
  PetscPrintf(PETSC_COMM_WORLD,
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

//__INSERT_LICENSE__
//$Id: compdate.cpp,v 1.1 2002/09/10 23:47:53 mstorti Exp $
#include <cstdio>
#include <petsc.h>

void print_petscfem_link_date() {
  PetscPrintf(PETSC_COMM_WORLD,
	      "NS-module built on " __DATE__ ", " __TIME__ ", with "
#ifdef __GNUC__
	      "GCC %d.%d.%d\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__
#else
	      "<unknown compiler>\n"
#endif
	      );
}



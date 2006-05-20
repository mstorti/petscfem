//__INSERT_LICENSE__
// $Id: generror.cpp,v 1.4.36.1 2006/05/20 21:11:19 dalcinl Exp $

#include <string>

using namespace std;

#include <petsc.h>
#include <src/generror.h>
#include <src/autostr.h>

GenericError PETSCFEM_GENERIC_ERROR;

GenericError::
GenericError(char *s,va_list ap) {
  AutoString as;
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}

GenericError::
GenericError(char *s,...) {
  AutoString as;
  va_list ap;

  va_start(ap,s);
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}

extern int MY_RANK,SIZE;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void petscfem_check_par_err(int ierro,GenericError &ge) {
  int ierr = MPI_Bcast (&ierro,1,MPI_INT,0,PETSCFEM_COMM_WORLD);	
#if 0
  PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			  "[%d] ierro %d\n",ierro,MY_RANK);
  PetscSynchronizedFlush(PETSCFEM_COMM_WORLD); 
  PetscFinalize();
  exit(0);
#endif
  if (ierro) {
    PetscPrintf(PETSCFEM_COMM_WORLD,"%s",ge.c_str());
    PetscFinalize();
    exit(0);
  }
}

//__INSERT_LICENSE__
// $Id: generror.cpp,v 1.5.10.1 2007/02/19 20:23:56 mstorti Exp $

#include <cstdlib>
#include <string>

using namespace std;

#include <petsc.h>
#include <src/fem.h>
#include <src/generror.h>
#include <src/autostr.h>

extern MPI_Comm PETSCFEM_COMM_WORLD;
extern int MY_RANK,SIZE;


GenericError PETSCFEM_GENERIC_ERROR;

GenericError::
GenericError(const char *s,va_list ap) {
  AutoString as;
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}

GenericError::
GenericError(const char *s,...) {
  AutoString as;
  va_list ap;

  va_start(ap,s);
  as.vsprintf(s,ap);
  *this = GenericError(string(as.str()));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void petscfem_check_par_err(int ierro,GenericError &ge) {
  int PFUNUSED ierr = MPI_Bcast (&ierro,1,MPI_INT,0,PETSCFEM_COMM_WORLD);	
#if 0
  ierr = PetscSynchronizedPrintf(PETSCFEM_COMM_WORLD,
			  "[%d] ierro %d\n",ierro,MY_RANK);
  ierr = PetscSynchronizedFlush(PETSCFEM_COMM_WORLD); 
  ierr = PetscFinalize();
  exit(0);
#endif
  if (ierro) {
    ierr = PetscPrintf(PETSCFEM_COMM_WORLD,"%s",ge.c_str());
    ierr = PetscFinalize();
    exit(0);
  }
}

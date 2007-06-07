// $Id$

#ifndef PF4PY_GVARS_H
#define PF4PY_GVARS_H

#include "Comm.h"
#include "Options.h"
#include "Mesh.h"
#include "Dofset.h"

extern TextHashTable* GLOBAL_OPTIONS;
extern Mesh*          GLOBAL_MESH;
extern int            SIZE, MY_RANK;

PF4PY_NAMESPACE_BEGIN

struct GlobalVars {
  GlobalVars(MPI_Comm       comm,
	     Options::Impl* thash  = NULL,
	     Mesh::Impl*    mesh   = NULL,
	     Dofset::Impl*  dofmap = NULL)
  {
    // set global variables used by PETScFEM
    PETSCFEM_COMM_WORLD = comm;
    MPI_Comm_size(PETSCFEM_COMM_WORLD, &SIZE);
    MPI_Comm_rank(PETSCFEM_COMM_WORLD, &MY_RANK);
    GLOBAL_MESH = mesh;
    if (thash)
      GLOBAL_OPTIONS = thash;
    else
      GLOBAL_OPTIONS = Options::GLOBALS;
    GLOBAL_OPTIONS->set_as_global();
    // try to clear any specific error message in PETSc
    char *msg;
    PetscErrorMessage(0, PETSC_NULL, &msg);
    if(msg) msg[0] = '\0';
  }
  ~GlobalVars()
  {
    // restore global variables used by PETScFEM
    PETSCFEM_COMM_WORLD = MPI_COMM_NULL;
    GLOBAL_MESH         = NULL;
    GLOBAL_OPTIONS      = Options::GLOBALS;
    GLOBAL_OPTIONS->set_as_global();
  }
};

PF4PY_NAMESPACE_END

#endif // PF4PY_GVARS_H

// Local Variables:
// mode: C++
// End:

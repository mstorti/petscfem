// $Id$

#ifndef PF4PY_MACROS_H
#define PF4PY_MACROS_H


#define PF4PY_INCREF(obj) do {if((obj)) (obj)->incref();} while(0)
#define PF4PY_DECREF(obj) do {if((obj)) (obj)->decref();} while(0)

#define PF4PY_ASSERT(cond, message) \
do {if (!(cond)) throw Error((message));} while(0)


#define  PF4PY_DELETE(delop, member) \
do {if ((member) != NULL) { delop ((member)); (member) = NULL; }} while(0)

#define  PF4PY_DELETE_FUNC(delop, member) \
do {if ((member) != NULL) { delop ((member)); (member) = NULL; }} while(0)

#define  PF4PY_DELETE_SCLR(member) PF4PY_DELETE_FUNC(delete,   member)
#define  PF4PY_DELETE_VCTR(member) PF4PY_DELETE_FUNC(delete[], member)

#define  PF4PY_PETSC_DESTROY(destroy, object)\
do {\
  if ((object) != PETSC_NULL && !PetscFinalizeCalled) {\
    PetscErrorCode ierr = destroy ((object)); (object) = NULL;\
    if (ierr) throw Error("PETSc error, calling function " #destroy);\
    }\
} while(0)

#define  PF4PY_PETSC_CALL(function, args)\
do {\
  PetscErrorCode ierr = function args;\
  if (ierr) throw Error(ierr, "PETSc error calling " #function);\
} while(0)


#endif // PF4PY_MACROS_H

// Local Variables:
// mode: C++
// End:

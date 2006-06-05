// $Id: macros.h,v 1.1.2.9 2006/06/05 16:00:54 dalcinl Exp $

#ifndef PYPF_MACROS_H
#define PYPF_MACROS_H


#define SMARTPTR(CLASS) public SmartPtr< ::CLASS >,

#define REFCOUNTER(CLASS) : public RefCounter

#define PYPF_INCREF(obj) if((obj)) (obj)->incref()
#define PYPF_DECREF(obj) if((obj)) (obj)->decref()

#define PYPF_ASSERT(cond, message) \
do {if (!(cond)) throw Error((message));} while(0)


#define  PYPF_DELETE(delop, member) \
do {if ((member) != NULL) { delop ((member)); (member) = NULL; }} while(0)

#define  PYPF_DELETE_FUNC(delop, member) \
do {if ((member) != NULL) { delop ((member)); (member) = NULL; }} while(0)

#define  PYPF_DELETE_SCLR(member) PYPF_DELETE(delete,   member)
#define  PYPF_DELETE_VCTR(member) PYPF_DELETE(delete[], member)

#define  PYPF_PETSC_DESTROY(destroy, object)\
do {\
  if ((object) != PETSC_NULL && !PetscFinalizeCalled) {\
    PetscErrorCode ierr = destroy ((object)); (object) = NULL;\
    if (ierr) throw Error("PETSc error, calling function " #destroy);\
    }\
} while(0)

#define  PYPF_PETSC_CALL(function, args)\
do {\
  PetscErrorCode ierr = function args;\
  if (ierr) throw Error("PETSc error, calling function " #function);\
} while(0)


#endif // PYPF_MACROS_H

// Local Variables:
// mode: C++
// End:

// $Id: macros.h,v 1.1.2.5 2006/04/27 19:09:17 rodrigop Exp $

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

#define  PYPF_PETSC_DESTROY(destroy, object) \
do { \
  if ((object) != PETSC_NULL && !PetscFinalizeCalled) \
    { destroy ((object)); (object) = NULL; } \
} while(0)

#endif // PYPF_MACROS_H

// Local Variables:
// mode: C++
// End:

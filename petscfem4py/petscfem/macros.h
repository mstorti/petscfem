// -*- c++ -*-
// $Id: macros.h,v 1.1.2.2 2006/03/20 16:06:00 rodrigop Exp $

#ifndef PYPF_MACROS_H
#define PYPF_MACROS_H


// #define PYPF_OBJ_GETOPTTBL_DECL \
// protected: \
// OptionTable* get_opt_table() const; \
// private:

#define SMARTPTR(CLASS) public SmartPtr< ::CLASS >,

#define  PYPF_DELETE(delop, member) \
do {if ((member) != NULL) { delop ((member)); (member) = NULL; }} while(0)

#define  PYPF_DELETE_SCLR(member) PYPF_DELETE(delete,   member)
#define  PYPF_DELETE_VCTR(member) PYPF_DELETE(delete[], member)

#define  PYPF_PETSC_DESTROY(destroy, object) \
do { \
  if ((object) != PETSC_NULL && !PetscFinalizeCalled) \
    { destroy ((object)); (object) = NULL; } \
} while(0)

#endif // PYPF_MACROS_H

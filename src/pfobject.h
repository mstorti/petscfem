// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: pfobject.h,v 1.5 2003/03/06 22:45:59 mstorti Exp $
#ifndef PETSCFEM_PFOBJECT_H
#define PETSCFEM_PFOBJECT_H

#include <string>
#include <src/fem.h>
#include <src/dofmap.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This is the most basic PETSc-FEM object. It reads itself from
    a stream (through a #FileStack#). */ 
class BasicObject {
public:
  /** The object reads itself from the #FileStack# and may modify
      #mesh# or #dofmap# through providing constraints, etc... 
      @param fstack (input) data is read from here
      @param mesh (input/output) as a side effect this structure may be modified
      @param dofmap (input/output) as a side effect this structure may be 
      modified. For instance it can add constraints or fixations. */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {}
  /** Dtor. */ 
  virtual ~BasicObject()=0;
  /** Produces objects, probably from a string identifying 
      the specific type. */ 
  static BasicObject *factory(string &type);
};

typedef BasicObject *BasicObject_factory_t(string &type);
extern BasicObject_factory_t *BasicObject_application_factory;

#endif

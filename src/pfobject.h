// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: pfobject.h,v 1.3 2003/02/25 20:34:22 mstorti Exp $
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

#endif

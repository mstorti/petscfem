// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: pfobject.h,v 1.2 2003/02/25 13:31:41 mstorti Exp $
#ifndef PETSCFEM_PFOBJECT_H
#define PETSCFEM_PFOBJECT_H

#include <string>
#include <src/fem.h>
#include <src/dofmap.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class BasicObject {
public:
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap) {}
  virtual ~BasicObject()=0;
  static BasicObject *factory(string &type);
};

#endif

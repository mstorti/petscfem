// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: pfobject.h,v 1.1 2003/02/25 04:18:10 mstorti Exp $
#ifndef PETSCFEM_PFOBJECT_H
#define PETSCFEM_PFOBJECT_H

#include <string>
#include <src/fstack.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class BasicObject {
public:
  virtual void read(FileStack *fstack) {}
  virtual ~BasicObject()=0;
  static BasicObject *factory(string &type);
};

#endif

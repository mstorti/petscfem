// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nullvort.h,v 1.1 2003/02/25 20:34:22 mstorti Exp $
#ifndef PETSCFEM_NULLVORT_H
#define PETSCFEM_NULLVORT_H

#include <src/pfobject.h>

/** This objects provides contranints that restrict the velocity
    field at a certain set of nodes to be irrotational.  */ 
class null_vort : public BasicObject {
private:
  /// The options table
  TextHashTable thash;
public:
  /** Dtor.  */ 
  ~null_vort() {}
  /** Specific read function for this obkect.  */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap);
};

#endif

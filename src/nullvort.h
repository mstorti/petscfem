// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nullvort.h,v 1.2 2003/02/26 01:46:33 mstorti Exp $
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
  /** Ctor.  */ 
  null_vort();
  /** Dtor.  */ 
  ~null_vort();
  /** Specific read function for this obkect.  */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap);
};

#endif

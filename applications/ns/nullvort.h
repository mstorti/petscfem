// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nullvort.h,v 1.1 2003/03/06 20:49:04 mstorti Exp $
#ifndef PETSCFEM_NULLVORT_H
#define PETSCFEM_NULLVORT_H

#include <src/pfobject.h>

/** This objects provides contranints that restrict the velocity
    field at a certain set of nodes to be irrotational.  */ 
class null_vort_bo : public BasicObject {
private:
  /// The options table
  TextHashTable thash;
  
public:
  /** Ctor.  */ 
  null_vort_bo();
  /** Dtor.  */ 
  ~null_vort_bo();
  /** Specific read function for this obkect.  */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap);
};

#endif

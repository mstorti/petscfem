// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nullvort.h,v 1.2 2003/03/07 03:13:06 mstorti Exp $
#ifndef PETSCFEM_NULLVORT_H
#define PETSCFEM_NULLVORT_H

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/pfobject.h>
#include "./nsi_tet.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This objects provides contranints that restrict the velocity
    field at a certain set of nodes to be irrotational.  */ 
class null_vort_bo : public BasicObject {
public:
  /** Ctor.  */ 
  null_vort_bo();
  /** Dtor.  */ 
  ~null_vort_bo();
  /** Specific read function for this obkect.  */ 
  virtual void read(FileStack *fstack,Mesh *mesh,Dofmap *dofmap);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class null_vort : public NonLinearRes {
public:
  /// Number of restrictions (=2)
  int nres() {return 2;};
  /// Initialize the elemset (maybe reads hash table)
  void init() {}
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
    FastMat2 & jac) {}
  /// Contructor
  null_vort() {};
  /// Destructor
  ~null_vort() {};
};

#endif

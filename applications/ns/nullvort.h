// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: nullvort.h,v 1.3 2003/03/07 21:23:52 mstorti Exp $
#ifndef PETSCFEM_NULLVORT_H
#define PETSCFEM_NULLVORT_H

#include <src/fem.h>
#include <src/readmesh.h>
#include <src/pfobject.h>
#include <src/cloud2.h>
#include <applications/ns/nsi_tet.h>
#include <applications/ns/lagmul.h>

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
/// Restriction for linearized free surface boundary condition
class null_vort : public LagrangeMult {
private:
  int fic_dof, nx, ndim;
  Cloud2 cloud;
  // coorinates of nodes in stencil, coordinates
  // of node, computed coefs.
  FastMat2 xlocc, x, x0, w, ww;
public:
  friend class null_vort_bo;
  /// Number of restrictions
  int nres() {return 1;};
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. */
  void lag_mul_dof(int jr,int &node,int &dof);
  /// Initialize the elemset (maybe reads hash table)
  void init();
  /// computes the residual and jacobian of the function to be imposed
  void res(int k,FastMat2 &U,FastMat2 & r,FastMat2 & lambda,
	   FastMat2 & jac);
};

#if 0
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

#endif

// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: advabso.h,v 1.7 2005/01/28 18:16:44 mstorti Exp $
#ifndef PETSCFEM_ADVABSO_H
#define PETSCFEM_ADVABSO_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#include "./lagmul.h"
#include "./gasflow.h"

#define gasflow_abso gasflow_abso2

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class AdvectiveAbso : public AdvdifLagrangeMult {
private:
  int ndim,nel,ndof;
  double Rgas;
  // Pointer to adv-diff flux fun
  NewAdvDifFF *adv_diff_ff;
  FastMat2 dummy,flux,fluxd,A_grad_U,
    grad_U,Uold,normal,A_jac,S,invS,c,
    Pi_m,Pi_p,Uref,tmp1,Ulambda,Uo,
    dU,Cp,invCp;
  int use_old_state_as_ref;
  // Per node normal
  Property normal_prop;
public:
  AdvectiveAbso(NewAdvDifFF *ff) 
    : adv_diff_ff(ff),
      use_old_state_as_ref(0) {} 
  ~AdvectiveAbso() { delete adv_diff_ff; } 
  int nres() { return ndof; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    node = 2; dof=jr;
  }
  //element init
  virtual void 
  element_hook(ElementIterator &element);
  void lm_initialize();
  void init();
  void res(int k,FastMat2 &U,FastMat2 &r,
	   FastMat2 &w,FastMat2 &jac);
  void close() {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class gasflow_abso : public AdvectiveAbso {
public:
  gasflow_abso() 
    :  AdvectiveAbso(new gasflow_ff(this)) { } 
};

#undef gasflow_abso

#endif

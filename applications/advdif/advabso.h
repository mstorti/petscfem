// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: advabso.h,v 1.1 2005/01/25 23:27:53 mstorti Exp $
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
  int ndim;
  double Rgas;
  // Poiter to adv-diff flux fun
  NewAdvDifFF *adv_diff_ff;
  FastMat2 dummy,flux,fluxd,A_grad_U,
    grad_U;
public:
  AdvectiveAbso(NewAdvDifFF *ff) 
    : adv_diff_ff(ff) {} 
  ~AdvectiveAbso() { delete adv_diff_ff; } 
  int nres() { return ndof; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    node = 2; dof=jr;
  }
  void lm_initialize() { }
  void init() {
    int ierr;
    TGETOPTDEF_ND(thash,int,ndim,0);
    flux.resize(1,ndof);
    fluxd.resize(1,ndof);
    A_grad_U.resize(1,ndof);
    grad_U.resize(2,ndim,ndof).set(0.);
  }
  void res(int k,FastMat2 &U,FastMat2 &r,
	   FastMat2 &w,FastMat2 &jac) {
    double delta_sc=0.0,
      lambda_max_pg=0.0;
    adv_diff_ff
      ->compute_flux(U, dummy, dummy, dummy, flux, fluxd,
		     A_grad_U, grad_U, dummy,
		     dummy, delta_sc, lambda_max_pg, dummy,
		     dummy, dummy, dummy, 0);
  }
  void close() {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class gasflow_abso : public AdvectiveAbso {
public:
  gasflow_abso() 
    :  AdvectiveAbso(new gasflow_ff()) { } 
};

#endif

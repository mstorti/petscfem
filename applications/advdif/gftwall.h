// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: gftwall.h,v 1.4 2005/01/27 00:48:08 mstorti Exp $
#ifndef PETSCFEM_GFTWALL_H
#define PETSCFEM_GFTWALL_H

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <src/lagmul.h>
#include "./lagmul.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Impose Dirichlet b.c.'s via Lagrange Multipliers. */ 
class gasflow_twall : public AdvdifLagrangeMult {
private:
  int ndim;
  double Rgas;
public:
  ~gasflow_twall() {} 
  int nres() { return 1; }
  void lag_mul_dof(int jr,int &node,int &dof) {
    node = 2; dof=1;
  }
  void lm_initialize() { }
  void init() {
    int ierr;
    NSGETOPTDEF_ND(int,ndim,0);
    //o Constant of a particular gas for ideal gas
    //  law (state equation for the gas)
    NSGETOPTDEF_ND(double,Rgas,287.);
  }
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac) {
    double rho = U.get(1,1);
    double p = U.get(1,ndim+2);
    double T = U.get(2,2);
    r.setel(p-rho*Rgas*T,1);
    w.set(0.);
    w.setel(1.0,1,ndim+2,1);
    jac.setel(-Rgas*T,1,1,1);
    jac.setel(1.0,1,1,ndim+2);
  }
  void close() {}
};

#endif

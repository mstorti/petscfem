// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: advec.h,v 1.1 2002/07/11 00:34:32 mstorti Exp $
#ifndef ADVEC_H
#define ADVEC_H

#include "advective.h"

class advec_ff : public NewAdvDifFF {
  int ndim;
  double diff;
  FastMat2 u,A,D,uu,tmp1;
 public:
  void start_chunk(int &options) {
    options &= ~SCALAR_TAU;
    ndim = 2;
    u.resize(1,ndim);
    u.set(0.).setel(1.,1);
    uu = sqrt(u.sum_square_all());
    A.resize(1,ndim).set(u).reshape(3,ndim,1,1);
    D.resize(2,ndim,ndim).eye(diff).reshape(4,ndim,ndim,1,1);
  }
  void element_hook(ElementIterator &element) {}
  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
    A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
  }
  virtual void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    // Should be used only for Absorbing boundary conditions
    assert(0);
  }
  void compute_flux(const FastMat2 &U,
		    const FastMat2 &iJaco, FastMat2 &H,
		    FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
		    FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
		    FastMat2 &tau_supg, double &delta_sc,
		    double &lam_max,FastMat2 &nor, FastMat2 &lambda,
		    FastMat2 &Vr, FastMat2 &Vr_inv,int options) {
    flux.prod(A,U,2,1,-1,-1);
    fluxd.prod(D,grad_U,1,-1,2,-2,-1,-2);
    A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
    G_source.set(0.);
    tau_supg.set(0.);
    delta_sc = 0.;
    lam_max = uu;
  }
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & grad_N,double w) {
    tmp1.set(dshapex).scale(w*diff);
    grad_N_D_grad_N.ir(2,1).ir(4,1);
    grad_N_D_grad_N.prod(tmp1,dshapex,-1,1,-1,2);
    grad_N_D_grad_N.rs();
  }
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    N_N_C.set(0.);
  }
  comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	     FastMat2 &N,double w) {
    N_P_C.set(0.);
  }
  void comp_P_supg(FastMat2 &P_supg) {
    tmp2.prod(A,grad_N(),-1,2,3,-1,1);
    P_supg.prod(tmp2,,1,-1,-1,2);
  }
  
};


#endif

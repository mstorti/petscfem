// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: advec.h,v 1.3 2002/07/11 16:35:53 mstorti Exp $
#ifndef ADVEC_H
#define ADVEC_H

#include "advective.h"

class advec_ff : public NewAdvDifFF {
  int ndim;
  double diff,uu,tau,tau_fac;
  FastMat2 u,A,D,tmp1,Uintri,tmp0;
  const NewAdvDif* e;
 public:
  advec_ff(const NewAdvDif *e) : NewAdvDifFF(e) {}
  ~advec_ff() {}
  void start_chunk(int &options) {
    int ierr;
    options &= ~SCALAR_TAU;
    ndim = 2;
    u.resize(1,ndim);
    u.set(0.).setel(1.,1);
    diff = 1e-3;
    uu = sqrt(u.sum_square_all());
    A.resize(1,ndim).set(u).reshape(3,ndim,1,1);
    D.resize(2,ndim,ndim).eye(diff).reshape(4,ndim,ndim,1,1);
    e = dynamic_cast<const NewAdvDif *>(elemset); 
    Uintri.resize(1,ndim);
    EGETOPTDEF_ND(elemset,double,tau_fac,1.); //nd
    int nel,ndof,nelprops;
    elemset->elem_params(nel,ndof,nelprops);
    enthalpy_fun = &identity_ef;
    identity_ef.init(ndof,ndim,nel);
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
    fluxd.prod(D,grad_U,2,-1,1,-2,-1,-2);
    A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
    G_source.set(0.);
    tau_supg.set(0.);
    delta_sc = 0.;
    lam_max = uu;
    Uintri.prod(iJaco,u,1,-1,-1);
    double Uh = sqrt(Uintri.sum_square_all());
    FastMat2::branch();
    if (uu*uu > 1e-5*Uh*diff) { // remove singularity when v=0
      FastMat2::choose(0);
      double Pe  = uu*uu/(Uh*diff);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = tau_fac/Uh*magic; // intrinsic time
    } else {
      FastMat2::choose(1);
      double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
      tau = tau_fac*h*h/(12.*diff);
    }
    FastMat2::leave();
    tau_supg.setel(tau,1,1);
  }
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & grad_N,double w) {
    tmp1.set(grad_N).scale(w*diff);
    grad_N_D_grad_N.ir(2,1).ir(4,1);
    grad_N_D_grad_N.prod(tmp1,grad_N,-1,1,-1,2);
    grad_N_D_grad_N.rs();
  }
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    N_N_C.set(0.);
  }
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	     FastMat2 &N,double w) {
    N_P_C.set(0.);
  }
#if 0
  void comp_P_supg(FastMat2 &P_supg) {
    P_supg.prod(A,*e->grad_N(),-1,2,3,-1,1).scale(tau);
  }
#endif
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// The `advec' (river or channel) element.
class advec : public NewAdvDif {
public:
  advec() :  NewAdvDif(new advec_ff(this)) {};
};

#endif

// -*- mode: C++ -*- 
//__INSERT_LICENSE__
//$Id: burgers.h,v 1.1 2001/04/01 21:47:38 mstorti Exp $
#ifndef BURGERS_H
#define BURGERS_H

#include "nwadvdifj.h"

#if 0
class burgers_ff : public NewAdvDifFF {
private:  
  double diffusivity, tau_fac, flux_law_coefficient;
  FastMat2 u0,u,Uintri,tmp0;
  NullSourceTerm source_term_;
  NullCJac c_jac_;
  UGlobalAJac a_jac_;
  GlobalScalarDJac d_jac_;
  
public:
  burgers_ff(NewElemset *elemset_) : a_jac(&a_jac_), source_term(&source_term_), 
    c_jac(&c_jac_), d_jac(&d_jac_), elemset(elemset_) {};
  void start_chunk(int &ret_options);
  void element_hook(ElementIterator &element);
  void compute_flux(COMPUTE_FLUX_ARGS);
  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    a_jac->comp_A_jac_n(A_jac_n,normal);
  }
  void comp_A_grad_N(FastMat2 & A,FastMat2 & B) {
    a_jac->comp_A_grad_N(A,B);
  }
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 &dshapex,double w) {
    d_jac->comp_grad_N_D_grad_N(grad_N_D_grad_N,dshapex,w);
  }
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    c_jac->comp_N_N_C(N_N_C, N, w);
  }
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w) {
    c_jac->comp_N_P_C(N_P_C,P_supg,N,w);
  }
};

class newadvdif_burgers : public NewAdvDif {
public:
  newadvdif_burgers() : NewAdvDif(new burgers_ff(this)) {};
  // newadvdif_burgers() {adv_diff_ff = new newburgers_ff_t(this);};
};

class newbcconv_burgers : public NewBcconv {
public:
  newbcconv_burgers() : NewBcconv(new burgers_ff(this)) {};
};
#endif

class burgers_ff : public advdif_wjac_ff {
  double diffusivity, tau_fac, flux_law_coefficient;
  FastMat2 u0,u,Uintri,tmp0;
  int ndim,ndof,nel,nelprops;
 public:
  burgers_ff(NewElemset *elemset_) : advdif_wjac_ff(elemset_) {};
  void start_chunk(int &ret_options);
  void element_hook(ElementIterator &element) {};
  void compute_flux(COMPUTE_FLUX_ARGS);
};

class newadvdif_burgers : public NewAdvDif {
public:
  newadvdif_burgers() : NewAdvDif(new burgers_ff(this)) {};
};

class newbcconv_burgers : public NewBcconv {
public:
  newbcconv_burgers() : NewBcconv(new burgers_ff(this)) {};
};

#endif

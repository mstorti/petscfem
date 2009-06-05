// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: burgers.h,v 1.8 2005/02/23 01:40:34 mstorti Exp $
#ifndef PETSCFEM_BURGERS_H
#define PETSCFEM_BURGERS_H

#include "nwadvdifj.h"

class burgers_ff : public advdif_wjac_ff {
  double diffusivity, tau_fac, flux_law_coefficient;
  FastMat2 u0,Uintri,tmp0,Ucpy,iJaco_cpy;
  int ndim,ndof,nel,nelprops;
  UGlobalAJac u_a_jac;
  NullSourceTerm st;
  NullCJac n_c_jac;
  GlobalScalarDJac g_d_jac;
  IdentityEF identity_ef;
 public:
  burgers_ff(NewElemset *elemset_) 
    : advdif_wjac_ff(elemset_,&g_d_jac,&u_a_jac,&n_c_jac,&st) {};
  void start_chunk(int &ret_options);
  void set_state(const FastMat2 &U,const FastMat2 &grad_U) {}
  void element_hook(ElementIterator &element) {};
  void compute_flux(COMPUTE_FLUX_ARGS);
  void get_Ajac(FastMat2 &Ajac);
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

// -*- mode: C++ -*- 
#ifndef NWADVDIF_H
#define NWADVDIF_H

#include "advective.h"

/** The class AdvDif is a NewElemset class plus a
    advdif flux function object.
*/
class NewAdvDif : public NewElemset { 
  NewAdvDifFF *adv_diff_ff;
public:
  int ndim,ndof,nel;
  NewAdvDif(NewAdvDifFF *adv_diff_ff_) : adv_diff_ff(adv_diff_ff_) {};
  NewAssembleFunction new_assemble;
  ASK_FUNCTION;
};

class newadvecfm2_ff_t : public NewAdvDifFF {
private:  
  int shock_capturing,na,nd,nc;
  FastMat2 D_jac_l, C_jac_l, tmp0;
  double tau_fac;
  FastMat2 u,u2,Uintri,AA,Ucpy,iJaco_cpy,
    tmp2;
  vector<double> djacv,cjacv;
  double *djacvp,*cjacvp;
  ElementIterator element;
  Property advective_jacobians_prop;
  const double *advjac;
  int ndim,ndof,nel;
public:

  /// One velocity for all the fields
  class UGlobal;
  friend class UGlobal;
  class UGlobal : public AJac {
    FastMat2 tmp,tmp3;
    newadvecfm2_ff_t &ff;
  public:
    UGlobal(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_flux,comp_A_grad_U,comp_A_grad_N,
      comp_Uintri;
  };
  UGlobal u_global;

  /// One velocity per field
  class UPerField;
  friend class UPerField;
  class UPerField : public AJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp;
  public:
    UPerField(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_flux,comp_A_grad_U,comp_A_grad_N,
      comp_Uintri;
  };
  UPerField u_per_field;

  newadvecfm2_ff_t(NewAdvDif *elemset_) : NewAdvDifFF(elemset_),
    u_per_field(*this), u_global(*this) {};
  void start_chunk(int ret_options);
  void element_hook(ElementIterator &element);
  void compute_flux(COMPUTE_FLUX_ARGS);
};

class newadvdif_advecfm2 : public NewAdvDif {
public:
  newadvdif_advecfm2() : NewAdvDif(new newadvecfm2_ff_t(this)) {};
  // newadvdif_advecfm2() {adv_diff_ff = new newadvecfm2_ff_t(this);};
};

#endif

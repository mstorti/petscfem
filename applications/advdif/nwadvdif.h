// -*- mode: C++ -*- 
#ifndef NWADVDIF_H
#define NWADVDIF_H

#include "advective.h"

class AJac {
public:
  virtual void comp_A_grad_N(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_A_grad_U(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_Uintri(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_flux(FastMat2 & A,FastMat2 & B) =0 ;
};

class DJac {
public:
  virtual void comp_D_grad_N(FastMat2 & A,FastMat2 & B) =0 ;
  virtual void comp_fluxd(FastMat2 & A,FastMat2 & B) =0 ;
  virtual void comp_dif_per_field(FastMat2 &dif_per_field)=0;
};

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
  FastMat2 C_jac_l, tmp0;
  double tau_fac;
  FastMat2 u,u2,Uintri,AA,Ucpy,iJaco_cpy,
    tmp2,D_jac,dif_per_field;
  vector<double> djacv,cjacv;
  double *djacvp,*cjacvp;
  ElementIterator element;
  Property advective_jacobians_prop, diffusive_jacobians_prop;
  const double *advjac,*difjac;
  int ndim,ndof,nel;
  AJac *a_jac;
  DJac *d_jac;
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

  /// Full advective jacobian
  class FullAdvJac;
  friend class FullAdvJac;
  class FullAdvJac : public AJac {
    newadvecfm2_ff_t &ff;
  public:
    FullAdvJac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_flux,comp_A_grad_U,comp_A_grad_N,
      comp_Uintri;
  };
  FullAdvJac full_adv_jac;

  /// Full diffusive jacobian
  class FullDifJac;
  friend class FullDifJac;
  class FullDifJac : public DJac {
    newadvecfm2_ff_t &ff;
  public:
    FullDifJac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd,comp_D_grad_N;
    void comp_dif_per_field(FastMat2 &dif_per_field);
  };
  FullDifJac full_dif_jac;

  newadvecfm2_ff_t(NewAdvDif *elemset);
  void start_chunk(int ret_options);
  void element_hook(ElementIterator &element);
  void compute_flux(COMPUTE_FLUX_ARGS);
  void comp_A_grad_N(FastMat2 & A,FastMat2 & B) {
    a_jac->comp_A_grad_N(FastMat2 & A,FastMat2 & B);
  }
  void comp_D_grad_N(FastMat2 & A,FastMat2 & B) {
    d_jac->comp_D_grad_N(FastMat2 & A,FastMat2 & B);
  }
};

class newadvdif_advecfm2 : public NewAdvDif {
public:
  newadvdif_advecfm2() : NewAdvDif(new newadvecfm2_ff_t(this)) {};
  // newadvdif_advecfm2() {adv_diff_ff = new newadvecfm2_ff_t(this);};
};

#endif

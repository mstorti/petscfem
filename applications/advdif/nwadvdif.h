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

class AJac {
public:
  virtual void comp_A_grad_N(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_A_grad_U(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_Uintri(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_flux(FastMat2 & A,FastMat2 & B) =0 ;
};

class DJac {
public:
  virtual void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & dshapex,double w) =0 ;
  virtual void comp_fluxd(FastMat2 & A,FastMat2 & B) =0 ;
  virtual void comp_dif_per_field(FastMat2 &dif_per_field)=0;
  virtual void update(const double *difjac) {};
};

class newadvecfm2_ff_t : public NewAdvDifFF {
private:  
  int shock_capturing,na,nd,nc;
  FastMat2 C_jac_l, tmp0;
  double tau_fac;
  FastMat2 u,u2,Uintri,AA,Ucpy,iJaco_cpy,
    tmp2,D_jac,dif_per_field,tmp3,eye_ndof;
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
    FastMat2 D_grad_N;
  public:
    FullDifJac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd;
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void update(const double *difjac) {ff.D_jac.set(difjac);}
    void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			      FastMat2 & dshapex,double w);
  };
  FullDifJac full_dif_jac;

  /// A global tensor (the same for all fields)
  class GlobalDifTensor;
  friend class GlobalDifTensor;
  class GlobalDifTensor : public DJac {
    newadvecfm2_ff_t &ff;
    FastMat2 D_grad_N;
    FastMat2 tmp1,tmp2;
  public:
    GlobalDifTensor(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd;
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void update(const double *difjac) {ff.D_jac.set(difjac);}
    void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			      FastMat2 & dshapex,double w);
  };
  GlobalDifTensor global_dif_tensor;

  /// A global tensor (the same for all fields)
  class PerFieldDifTensor;
  friend class PerFieldDifTensor;
  class PerFieldDifTensor : public DJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp,tmp2;
  public:
    PerFieldDifTensor(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd;
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void update(const double *difjac) {ff.D_jac.set(difjac);}
    void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			      FastMat2 & dshapex,double w);
  };
  PerFieldDifTensor per_field_dif_tensor;

  /// Global scalar diffusion
  class GlobalScalar;
  friend class GlobalScalar;
  class GlobalScalar : public DJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp;
  public:
    GlobalScalar(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd;
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			      FastMat2 & dshapex,double w);
  };
  GlobalScalar global_scalar_djac;

  /// Scalar diffusion per field
  class ScalarDifPerField;
  friend class ScalarDifPerField;
  class ScalarDifPerField : public DJac {
    newadvecfm2_ff_t &ff;
    FastMat2 grad_N_grad_N,tmp;
  public:
    ScalarDifPerField(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    FastMat2Shell comp_fluxd;
    void update(const double *difjac) {ff.D_jac.set(difjac);}
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			      FastMat2 & dshapex,double w);
  };
  ScalarDifPerField scalar_dif_per_field;

  newadvecfm2_ff_t(NewAdvDif *elemset);
  void start_chunk(int ret_options);
  void element_hook(ElementIterator &element);
  void compute_flux(COMPUTE_FLUX_ARGS);
  void comp_A_grad_N(FastMat2 & A,FastMat2 & B) {
    a_jac->comp_A_grad_N(A,B);
  }
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 &dshapex,double w) {
    d_jac->comp_grad_N_D_grad_N(grad_N_D_grad_N,dshapex,w);
  }
};

class newadvdif_advecfm2 : public NewAdvDif {
public:
  newadvdif_advecfm2() : NewAdvDif(new newadvecfm2_ff_t(this)) {};
  // newadvdif_advecfm2() {adv_diff_ff = new newadvecfm2_ff_t(this);};
};

#endif

// -*- mode: C++ -*- 
#ifndef NWADVDIF_H
#define NWADVDIF_H

#include "advective.h"

class CJac {
public:
  virtual void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w)=0;
  virtual void comp_G_source(FastMat2 &G_source, FastMat2 &U)=0;
  virtual void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			  FastMat2 &N,double w)=0;
};

class AJac {
public:
  virtual void comp_A_grad_N(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_A_grad_U(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_Uintri(FastMat2 & A,FastMat2 & B)=0;
  virtual void comp_flux(FastMat2 & A,FastMat2 & B) =0 ;
  // virtual void update(const double *advjac) {ff.A_jac.set(advjac);}
};

class DJac {
public:
  virtual void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & dshapex,double w) =0 ;
  virtual void comp_fluxd(FastMat2 & A,FastMat2 & B) =0 ;
  virtual void comp_dif_per_field(FastMat2 &dif_per_field)=0;
  virtual void update(const double *difjac)=0;
};

class newadvecfm2_ff_t : public NewAdvDifFF {
private:  
  int shock_capturing,na,nd,nc;
  FastMat2 tmp0;
  double tau_fac;
  FastMat2 u,u2,Uintri,AA,Ucpy,iJaco_cpy,
    tmp2,D_jac,dif_per_field,tmp3,eye_ndof,
    C_jac,N_C;
  vector<double> djacv,cjacv;
  double *djacvp,*cjacvp;
  ElementIterator element;
  Property advective_jacobians_prop, 
    diffusive_jacobians_prop, reactive_jacobians_prop;
  const double *advjac,*difjac,*reacjac;
  int ndim,ndof,nel;
  AJac *a_jac;
  DJac *d_jac;
  CJac *c_jac;
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

  /// Full reactive Jacobian
  class FullCJac;
  friend class FullCJac;
  class FullCJac : public CJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp,tmp2,tmp26,tmp27;
  public:
    FullCJac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);
    void comp_G_source(FastMat2 &G_source, FastMat2 &U);
    void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		    FastMat2 &N,double w);
  };
  FullCJac full_c_jac;

  /// Global scalar reactive Jacobian
  class ScalarCJac;
  friend class ScalarCJac;
  class ScalarCJac : public CJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp,tmp2,tmp26,tmp27;
  public:
    ScalarCJac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);
    void comp_G_source(FastMat2 &G_source, FastMat2 &U);
    void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		    FastMat2 &N,double w);
  };
  ScalarCJac scalar_c_jac;

  /// Global scalar reactive Jacobian
  class ScalarPerFieldCjac;
  friend class ScalarPerFieldCjac;
  class ScalarPerFieldCjac : public CJac {
    newadvecfm2_ff_t &ff;
    FastMat2 tmp,tmp2,N_N_C,tmp26;
  public:
    ScalarPerFieldCjac(newadvecfm2_ff_t &ff_) : ff(ff_) {};
    void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);
    void comp_G_source(FastMat2 &G_source, FastMat2 &U);
    void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		    FastMat2 &N,double w);
  };
  ScalarPerFieldCjac scalar_per_field_c_jac;

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
    void update(const double *difjac) {};
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
    void comp_dif_per_field(FastMat2 &dif_per_field);
    void update(const double *difjac) {ff.D_jac.set(difjac);}
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
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    c_jac->comp_N_N_C(N_N_C, N, w);
  }
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w) {
    c_jac->comp_N_P_C(N_P_C,P_supg,N,w);
  }
};

class newadvdif_advecfm2 : public NewAdvDif {
public:
  newadvdif_advecfm2() : NewAdvDif(new newadvecfm2_ff_t(this)) {};
  // newadvdif_advecfm2() {adv_diff_ff = new newadvecfm2_ff_t(this);};
};

#endif

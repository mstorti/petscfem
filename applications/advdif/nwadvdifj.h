// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: nwadvdifj.h,v 1.2 2001/05/30 03:58:35 mstorti Exp $
#ifndef NWADVDIFJ_H
#define NWADVDIFJ_H

#include "nwadvdif.h"

/// Null source term
class NullSourceTerm : public SourceTerm {
 public:
  void add_source_term(FastMat2 &G_source) {};
};

class NullCJac : public CJac {
 public:
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    N_N_C.set(0.);}
  void comp_G_source(FastMat2 &G_source, FastMat2 &U) {
    G_source.set(0.);}
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w) {
    N_P_C.set(0.);}
};

class UGlobalAJac : public AJac {
  FastMat2 tmp,tmp2,tmp3,tmp5;
  int ndim,ndof,nel;
 public:
  FastMat2 u;
  FastMat2Shell comp_flux,comp_A_grad_U,comp_A_grad_N,
    comp_Uintri,comp_A_jac_n;
  void comp_vel_per_field(FastMat2 &vel_per_field);
  void init(int ndim_,int ndof_,int nel_);
};

class GlobalScalarDJac : public DJac {
  FastMat2 tmp,eye_ndof;
  int ndim,ndof,nel;
 public:
  double diff;
  FastMat2Shell comp_fluxd;
  void comp_dif_per_field(FastMat2 &dif_per_field);
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);
  void update(const double *difjac) {};
  void init(int ndim_,int ndof_,int nel_);
};

class advdif_wjac_ff : public NewAdvDifFF {
public:
  AJac *a_jac;
  DJac *d_jac;
  CJac *c_jac;
  SourceTerm *source_term;
  
  advdif_wjac_ff(NewElemset *elemset_,
		 DJac *d=NULL,AJac *a_jac_=NULL, 
		 CJac *c=NULL, SourceTerm *st=NULL);
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

#endif

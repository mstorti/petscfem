// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: nwadvdifj.h,v 1.5 2002/02/17 03:59:51 mstorti Exp $
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

/** In this class we have only to define the advective,
    diffusive and reactive jacobians, and source term.
*/
class advdif_wjac_ff : public NewAdvDifFF {

public:

  /// The advective jacobian
  AJac *a_jac;

  /// The diffusive jacobian
  DJac *d_jac;

  /// the reactive jacobian
  CJac *c_jac;

  /// The source term
  SourceTerm *source_term;

  /// Contructor from the objects
  advdif_wjac_ff(NewElemset *elemset_,
		 DJac *d=NULL,AJac *a_jac_=NULL, 
		 CJac *c=NULL, SourceTerm *st=NULL);

  /** Compute $!A\dot \nor = \sum_j A_j n_j$. Often, $\nor$ is a
      normal vector, but you can't assume that it is of unit length.
      @param A_jac_n (output) The contraction of the advective
      jacobian with the #nornal# vector.
      @param normal (input) The vector onto which to make the
      projection. 
   */ 
  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    a_jac->comp_A_jac_n(A_jac_n,normal);
  }

  /** Compute $\sum_{j} A_j \dep_j N$ where $N$ are the interpolation
      functions. 
      @param A_grad_N (output) The contraction.
      @param grad_N (input) the gradient of the interpolation function.
  */
  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
    a_jac->comp_A_grad_N(A_grad_N,grad_N);
  }

  /** Compute $\sum_{ij} w D_{ij} \dep_i N \dep_j N $. $w$ is a scalar weight,
      normally the weight of the Gauss point times the determinant of
      the transformation.
      @param grad_N_D_grad_N (output) the result.
      @param grad_N (input) the gradient of the interpolation
      functions
      @param w (input) the scalar weight
  */ 
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 &grad_N,double w) {
    d_jac->comp_grad_N_D_grad_N(grad_N_D_grad_N,grad_N,w);
  }
 
  /** Comput $\sum_{ij} C_{ij} N_i N_j$. 
      @param N_N_C (output) the result
      @param N (input) The interpolation function
      @param w (input) the scalar weight
  */ 
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    c_jac->comp_N_N_C(N_N_C, N, w);
  }
  
  /** Compute $\sum_{ij} C_{ij} N_i P_j$
      @param N_P_C (output) the result
      @param N (input) The interpolation function
      @param w (input) the scalar weight
  */ 
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w) {
    c_jac->comp_N_P_C(N_P_C,P_supg,N,w);
  }
};

#endif

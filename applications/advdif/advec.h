// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: advec.h,v 1.9 2002/12/14 21:19:01 mstorti Exp $
#ifndef ADVEC_H
#define ADVEC_H

#include "advective.h"

// This flag sets whether `P_supg' is computed in the flux function or 
// the standard computation P_supg = tau * A * grad_N is perfomed.
#define USE_COMP_P_SUPG

/** This a simple example of an advective class with
    computation of the weight function (more exactly of the
    perturbation function P_supg) not following the standard 
    rule "P_supg = tau * A * grad_N"
*/
class advec_ff : public NewAdvDifFF {
  /// dimension of the element
  int ndim;
  /// Diffusivity
  double diff;
  /// absolute value of velocity vector
  double uu;
  /// SUPG characteristic time
  double tau;
  /// Factor allowing user to scale the stabilization term
  double tau_fac;
  /// velocity vector (constant for the whole domain)
  FastMat2 u;
  /// Advective Jacobian
  FastMat2 A;
  /// Diffusive Jacobian
  FastMat2 D;
  /// Temporary matrices
  FastMat2 tmp1,tmp0;
  /// Velocity vector in intrinsic (master element) coordinates
  FastMat2 Uintri;
  /// Velocity property 
  Property u_prop;
  // Stores element number for use at the Gauss point level
  ElementIterator element_m;
 public:
  /// Constructor
  advec_ff(const NewAdvDif *e) : NewAdvDifFF(e) {}
  /// Destructor
  ~advec_ff() {}
  /// This is executed {\em before} the element loop.
  void start_chunk(int &options);
  void element_hook(ElementIterator &element) {
    element_m = element; }
  // A_grad_N(iel,mu,nu) = A(j,mu,nu) * grad_N(j,iel)
  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
    A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
  }
  // Should be used only for Absorbing boundary conditions
  virtual void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    assert(0);
  }
  void compute_flux(const FastMat2 &U,
		    const FastMat2 &iJaco, FastMat2 &H,
		    FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
		    FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
		    FastMat2 &tau_supg, double &delta_sc,
		    double &lam_max,FastMat2 &nor, FastMat2 &lambda,
		    FastMat2 &Vr, FastMat2 &Vr_inv,int options);
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & grad_N,double w);
  // No source term
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
    N_N_C.set(0.);
  }
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
	     FastMat2 &N,double w) {
    N_P_C.set(0.);
  }
#ifdef USE_COMP_P_SUPG
  void comp_P_supg(FastMat2 &P_supg) {
    P_supg.prod(A,*new_adv_dif_elemset->grad_N(),-1,2,3,-1,1).scale(tau);
  }
#endif
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// The `advec' element.
class advec : public NewAdvDif {
public:
  advec() :  NewAdvDif(new advec_ff(this)) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class advec_bcconv : public NewBcconv {
public:
  advec_bcconv() : NewBcconv(NULL) {};
};

#endif

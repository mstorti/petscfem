// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
//$Id: advec.h,v 1.8 2002/12/14 14:42:49 mstorti Exp $
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
 public:
  /// Constructor
  advec_ff(const NewAdvDif *e) : NewAdvDifFF(e) {}
  /// Destructor
  ~advec_ff() {}
  /// This is executed {\em before} the element loop.
  void start_chunk(int &options) {
    new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset); 
    int ierr;
    // Tell `advdife' that we will use a scalar `tau'
    options &= ~SCALAR_TAU;
    // Only 2D
    ndim = 2;
    // Resize velocity vector accordingly
    u.resize(1,ndim);
    u.set(0.).setel(1.,1);
    // read from options table
    ierr = elemset->get_double("u",*u.storage_begin(),1,ndim);
    // Set values
    EGETOPTDEF_ND(elemset,double,diff,1.);
    // Factor set by user scaling stabilization term
    EGETOPTDEF_ND(elemset,double,tau_fac,1.);
    // Absolute value of velocity
    uu = sqrt(u.sum_square_all());
    // Advective Jacobian
    A.resize(1,ndim).set(u).reshape(3,ndim,1,1);
    // Diffusive Jacobian
    D.resize(2,ndim,ndim).eye(diff).reshape(4,ndim,ndim,1,1);
    // Intrinsic velocity
    Uintri.resize(1,ndim);
    // Get element integer props
    int nel,ndof,nelprops;
    elemset->elem_params(nel,ndof,nelprops);
    // Set pointer to enthalpy function
    // Use identity entalphy for simplicity
    enthalpy_fun = &identity_ef;
    identity_ef.init(ndof,ndim,nel);
  }
  void element_hook(ElementIterator &element) {}
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
		    FastMat2 &Vr, FastMat2 &Vr_inv,int options) {
    // Convective flux
    // flux(j,mu) = A(j,mu,nu) * U(nu)
    flux.prod(A,U,2,1,-1,-1);
    // Diffusive flux
    // fluxd(j,mu) = D(j,k,mu,nu) * grad_U(k,nu)
    fluxd.prod(D,grad_U,2,-1,1,-2,-1,-2);
    // A_grad_U(mu) = A(j,mu,nu) * grad_U(j,nu)
    A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
    // Null source term 
    G_source.set(0.);
    // Set to null
    tau_supg.set(0.);
    // No shock capturing
    delta_sc = 0.;
    // maximum eigenvlue = absolute value of velocity
    lam_max = uu;
    // Intrinsic velocity
    Uintri.prod(iJaco,u,1,-1,-1);
    // This has scale of U/h
    double Uh = sqrt(Uintri.sum_square_all());
    FastMat2::branch();
    if (uu*uu > 1e-5*Uh*diff) {
      // remove singularity when v=0
      FastMat2::choose(0);
      double Pe  = uu*uu/(Uh*diff);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = tau_fac/Uh*magic; // intrinsic time
    } else {
      FastMat2::choose(1);
      // This id computed for very low velocity (for instance |u|=0)
      double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
      tau = tau_fac*h*h/(12.*diff);
    }
    FastMat2::leave();
    // Set tau_(1,1) = scalar tau
    tau_supg.setel(tau,1,1);
  }
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & grad_N,double w) {
    // tmp1 = w * diff * grad_N
    tmp1.set(grad_N).scale(w*diff);
    grad_N_D_grad_N.ir(2,1).ir(4,1);
    grad_N_D_grad_N.prod(tmp1,grad_N,-1,1,-1,2);
    grad_N_D_grad_N.rs();
  }
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

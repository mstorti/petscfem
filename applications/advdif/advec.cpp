//__INSERT_LICENSE__
//$Id: advec.cpp,v 1.3 2002/12/14 21:19:01 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "advec.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void advec_ff::start_chunk(int &options) {
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset); 
  int ierr;
  // Tell `advdife' that we will use a scalar `tau'
  options &= ~SCALAR_TAU;
  // Only 2D
  ndim = 2;
  // Resize velocity vector accordingly
  u.resize(1,ndim);
  u.set(0.).setel(1.,1);
  // get velocity property
  elemset->get_prop(u_prop,"u");
  assert(u_prop.length==ndim);
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

void advec_ff::compute_flux(const FastMat2 &U,
			    const FastMat2 &iJaco, FastMat2 &H,
			    FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
			    FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
			    FastMat2 &tau_supg, double &delta_sc,
			    double &lam_max,FastMat2 &nor, FastMat2 &lambda,
			    FastMat2 &Vr, FastMat2 &Vr_inv,int options) {
  // Set velocity vector
  u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void advec_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & grad_N,double w) {
  // tmp1 = w * diff * grad_N
  tmp1.set(grad_N).scale(w*diff);
  grad_N_D_grad_N.ir(2,1).ir(4,1);
  grad_N_D_grad_N.prod(tmp1,grad_N,-1,1,-1,2);
  grad_N_D_grad_N.rs();
}

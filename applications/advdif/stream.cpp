//__INSERT_LICENSE__

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "stream.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::start_chunk(int &ret_options) {
  v = 10.;
  int nel,ndof,nelprops;
  elemset->elem_params(nel,ndof,nelprops);
  IdentityEF *ef = dynamic_cast<IdentityEF *>(enthalpy_fun);
  assert(ef);
  ef->init(ndof,1,nel);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff
::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  options |= SCALAR_TAU;	// tell the advective element routine
				// that we are returning a scalar tau
  flux.ir(1,1).set(U).scale(v).rs();
  fluxd.set(0.);
  grad_U.ir(2,1);
  A_grad_U.set(grad_U).scale(v);
  grad_U.rs();

  if ( options & COMP_UPWIND ) {
    double h_supg = 2./iJaco.get(1,1);
    tau_supg.setel(h_supg/(2.*v),1,1);
    lam_max = v;
  }

  if (options & COMP_SOURCE)  G_source.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0); // Not implemented yet
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  grad_N.ir(1,1);
  A_grad_N.ir(2,1).ir(3,1).set(grad_N).scale(v);
  grad_N.rs();
  A_grad_N.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			  FastMat2 &dshapex,double w) {
  grad_N_D_grad_N.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		FastMat2 &N,double w) {
  N_P_C.set(0.);
}


//__INSERT_LICENSE__
// $Id: smoke.cpp,v 1.1 2003/05/25 13:52:05 mstorti Exp $

#include "./smoke.h"

void smoke_ff::start_chunk(int &ret_options) {
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);     

  // Get element integer props
  elemset->elem_params(nel,ndof,nelprops);

  elemset->get_prop(u_prop,"u");
  assert(u_prop.length==ndim);
  u.resize(1,ndim);
  // Tell `advdife' that we will use a scalar `tau'
  ret_options |= SCALAR_TAU;
  Cp.resize(2,ndof,ndof).eye();
  W_N.resize(2,nel,nel);
  A.resize(3,ndim,ndof,ndof);
}

void smoke_ff::element_hook(ElementIterator &element) {
  u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
  A.ir(2,1).ir(3,1).set(u).rs();
  element_m = element;
}

void smoke_ff::set_state(const FastMat2 &UU) { 
  phi = U.get(1); 
  U.set(UU);
}

void smoke_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0);			// fixme:= Not implemented yet, used for
				// absorbing boundary conditions
}

void smoke_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(A,grad_N,-1,2,3,-1,1);
}

void smoke_ff::compute_flux(COMPUTE_FLUX_ARGS) {
  u.set(new_adv_dif_elemset->prop_array(element_m,u_prop));
  double uu = sqrt(u.sum_square_all());
  // Convective flux
  // flux(j,mu) = A(j,mu,nu) * U(nu)
  flux.prod(A,U,2,1,-1,-1);
  // Diffusive flux
  // fluxd(j,mu) = D(j,k,mu,nu) * grad_U(k,nu)
  fluxd.set(0.);
  // A_grad_U(mu) = A(j,mu,nu) * grad_U(j,nu)
  A_grad_U.prod(A,grad_U,-1,1,-2,-1,-2);
  // Null source term (right now);
  G_source.set(0.);
  // Set to zero
  tau_supg.set(0.);
  // No shock capturing
  delta_sc = 0.;
  // maximum eigenvlue = absolute value of velocity
  lam_max = uu;
  // Intrinsic velocity
  Uintri.prod(iJaco,u,1,-1,-1);
  // This has scale of U/h
  double Uh = sqrt(Uintri.sum_square_all());

  // Consider pure advective flow
  double magic = 1., tau_fac=1.;

  // Set tau_(1,1) = scalar tau
  tau_supg.setel(tau_fac/Uh*magic,1,1);
}

void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			  FastMat2 & dshapex,double w) {
  grad_N_D_grad_N.set(0.);
}

void smoke_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.set(0.);
}

void smoke_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			  FastMat2 &N,double w) {
  N_P_C.set(0.);
}

void smoke_ff::enthalpy(FastMat2 &H) {  H.set(U); }

void smoke_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		 double w) {
  W_N.prod(W,N,1,2).scale(w);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

void smoke_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}


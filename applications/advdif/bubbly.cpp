//__INSERT_LICENSE__
//$Id: bubbly.cpp,v 1.1 2002/02/15 12:47:37 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::start_chunk(int &ret_options) {
  int ierr;
  elemset->elem_params(nel,ndof,nelprops);  
  EGETOPTDEF_ND(elemset,int,ndim,0);
  EGETOPTDEF_ND(elemset,double,rho_l,0);
  EGETOPTDEF_ND(elemset,double,rho_g,0);
  assert(ndim>0);
  assert(ndof==4+2*ndim);
  assert(rho_l>0.);
  assert(rho_g>0.);
  vl_indx = 3;
  vg_indx = vl_indx+ndim;
  k_indx = vg_indx+ndim;
  e_indx = k_indx+1;
  Cp.resize(2,ndof,ndof);
  Ajac.resize(3,ndim,ndof,ndof);
  Id.resize(2,ndim,ndim).eye();
  Amoml.resize(2,ndim,ndim);
  Amomg.resize(2,ndim,ndim);
  Y.resize(3,ndim,ndim,ndim);

  Djac.resize(4,ndim,ndof,ndim,ndof);
  tmp1.resize(4,ndim,nel,ndof,ndof);
  Cjac.resize(2,ndof,ndof);
  tmp2.resize(3,ndof,nel,ndof);
  tmp3.resize(2,ndof,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bubbly_ff::bubbly_ff(const NewAdvDif *e) 
  : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bubbly_ff::~bubbly_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::set_state(const FastMat2 &UU) {
  U.set(UU);
  // Scalar variables
  alpha_l = U.get(1);
  alpha_g = 1.-alpha_l;
  arho_l = alpha_l*rho_l;
  arho_g = alpha_g*rho_g;

  p = U.get(2);
  k = U.get(k_indx);
  eps = U.get(e_indx);
  p = U.get(2);
  // Velocities
  U.is(1,vl_indx,vl_indx+ndim-1);
  v_l.set(U);
  U.rs().is(1,vg_indx,vg_indx+ndim-1);
  v_g.set(U);
  U.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::enthalpy(FastMat2 &H) {
  // Scalar values
  H.setel(arho_l,1);
  H.setel(arho_g,2);
  H.setel(arho_l*k,k_indx);
  H.setel(arho_l*eps,e_indx);
  // Vector values
  H.is(1,vl_indx,vl_indx+ndim-1).set(v_l).scale(arho_l);
  H.rs();
  H.is(1,vg_indx,vg_indx+ndim-1).set(v_g).scale(arho_g);
  H.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff
::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  double arho_l,arho_g;

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Enthalpy Jacobian Cp
  // First column
  Cp.set(0.);
  Cp.setel(rho_l,1,1);
  Cp.setel(-rho_g,2,1);
  Cp.is(1,vl_indx,vl_indx+ndim-1).ir(2,1).set(v_l).scale(rho_l);
  Cp.is(1).is(1,vg_indx,vg_indx+ndim-1).set(v_g).scale(-rho_g);
  Cp.rs();
  Cp.setel(rho_l*k,k_indx,1);
  Cp.setel(rho_l*eps,e_indx,1);
  // V_l column block
  Cp.is(1,vl_indx,vl_indx+ndim-1).is(2,vl_indx,vl_indx+ndim-1)
    .eye(arho_l);
  // V_g column block
  Cp.rs().is(1,vg_indx,vg_indx+ndim-1).is(2,vg_indx,vg_indx+ndim-1)
    .eye(arho_g);
  // Turbulence
  Cp.rs();
  Cp.setel(arho_l,k_indx,k_indx);
  Cp.setel(arho_l,e_indx,e_indx);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Advective fluxes
  flux.ir(1,1).set(v_l).scale(rho_l*alpha_l);
  flux.ir(1,2).set(v_g).scale(rho_g*alpha_g);

  flux.rs();
  Amoml.prod(v_l,v_l,1,2).scale(rho_l)
    .axpy(Id,p);
  flux.is(1,vl_indx,vl_indx+ndim-1).axpy(Amoml,alpha_l);

  flux.rs();
  Amomg.prod(v_g,v_g,1,2).scale(rho_g)
    .axpy(Id,p);
  flux.is(1,vg_indx,vg_indx+ndim-1).axpy(Amomg,alpha_g);

  // TUrb
  flux.rs();
  flux.ir(1,k_indx).set(v_l).scale(rho_l*alpha_l*k);
  flux.ir(1,e_indx).set(v_l).scale(rho_l*alpha_l*eps);
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Adjective Jacobians
  // First column
  Ajac.set(0.);
  Ajac.ir(2,1).ir(3,1).set(v_l).scale(rho_l);
  Ajac.ir(2,2).set(v_g).scale(-rho_g);
  Ajac.rs().is(2,vl_indx,vl_indx+ndim-1).ir(3,1).set(Amoml);
  Ajac.rs().is(2,vg_indx,vg_indx+ndim-1).ir(3,1).axpy(Amomg,-1.);
  Ajac.rs().ir(2,k_indx).ir(3,1).set(v_l).scale(rho_l*k);
  Ajac.rs().ir(2,e_indx).ir(3,1).set(v_l).scale(rho_l*eps);
  
  // Second column
  Ajac.rs().is(2,vl_indx,vl_indx+ndim-1).ir(3,2).axpy(Id,alpha_l);
  Ajac.rs().is(2,vg_indx,vg_indx+ndim-1).ir(3,2).axpy(Id,alpha_g);
  
  Ajac.rs().ir(2,1).is(3,vl_indx,vl_indx+ndim-1).axpy(Id,arho_l);
  Ajac.rs().ir(2,2).is(3,vg_indx,vg_indx+ndim-1).axpy(Id,arho_g);
  
  /// fixme:= Verify this!!!
  Y.prod(v_l,Id,1,2,3).scale(arho_l);
  Ajac.rs().is(2,vl_indx,vl_indx+ndim-1).is(3,vl_indx,vl_indx+ndim-1)
    .add(Y);
  Y.exc(2,3);
  Ajac.add(Y);

  Y.rs().prod(v_g,Id,1,2,3).scale(arho_g);
  Ajac.rs().is(2,vg_indx,vg_indx+ndim-1).is(3,vg_indx,vg_indx+ndim-1)
    .add(Y);
  Y.exc(2,3);
  Ajac.add(Y);
  Y.rs();
  
  Ajac.rs().ir(2,k_indx).is(3,vl_indx,vl_indx+ndim-1)
    .axpy(Id,arho_l*k);
  Ajac.rs().ir(2,e_indx).is(3,vl_indx,vl_indx+ndim-1)
    .axpy(Id,arho_l*eps);

  Ajac.rs().ir(2,k_indx).ir(3,k_indx).set(v_l).scale(arho_l);
  Ajac.rs().ir(2,e_indx).ir(3,e_indx).set(v_l).scale(alpha_g*rho_l);

  Cjac.set(0.);
  if (options & COMP_SOURCE)  G_source.set(0.);

  Djac.reshape(2,ndim*ndof,ndim*ndof).eye().reshape(4,ndim,ndof,ndim,ndof);
  fluxd.prod(Djac,grad_U,1,2,-1,-2,-1,-2);
  tau_supg.set(0.).eye(1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0); // fixme:= not defined yet (for use with
	     // absorbing b.c's and bcconv)
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(Ajac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				     FastMat2 &dshapex,double w) {
  tmp1.prod(Djac,dshapex,-1,2,3,4,-1,1).scale(w);
  grad_N_D_grad_N.prod(Djac,dshapex,1,2,-1,4,-1,3);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,Cjac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,Cjac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}


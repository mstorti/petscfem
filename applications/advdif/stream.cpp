//__INSERT_LICENSE__
//$Id: stream.cpp,v 1.15 2002/02/19 18:08:59 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "stream.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DummyEnthalpyFun
::set_state(const FastMat2 &U) { 
  s->set_state(U); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DummyEnthalpyFun
::enthalpy(FastMat2 &H) {
  s->enthalpy(H); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DummyEnthalpyFun
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double w) {
  s->comp_W_Cp_N(W_Cp_N,W,N,w);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DummyEnthalpyFun
::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  s->comp_P_Cp(P_Cp,P_supg);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::start_chunk(int &ret_options) {
  // Initialize function objects
  if (!channel) channel = ChannelShape::factory(elemset);
  channel->init();

  if (!friction_law) 
    friction_law = FrictionLaw::factory(elemset);
  friction_law->init();

  // elemset->get_prop(slope_prop,"slope");
  // assert(slope_prop.length==1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::element_hook(ElementIterator &element) {
  // S = elemset->prop_val(element,slope_prop);
  channel->element_hook(element);
  friction_law->element_hook(element);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stream_ff::stream_ff(const NewAdvDif *e) 
  : AdvDifFFWEnth(e), channel(NULL), 
  friction_law(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stream_ff::~stream_ff() { delete friction_law; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::set_state(const FastMat2 &U) {
  u = U.get(1);
  channel->geometry(u,area,wl_width,perimeter);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::enthalpy(FastMat2 &H) {
  H.setel(area,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_Cp_N.ir(2,1).ir(4,1);
  W_Cp_N.prod(W,N,1,2).scale(wl_width*weight).rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.set(P_supg).scale(wl_width);
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

  // The equation is dA/dt + dQ/dx = 0
  // 
  // A:= the transversal area is equivalent to the enthalpy (or the
  //               conservative variable)
  // Q:= the volumetric flow, is equivalent to the flux
  // 
  // Q and A are functions of the unknown `u(x)' where
  // u:= is the depth of the fluid
  // `Q=Q(A)' through the `friction law', and `A=A(u)' through the
  // geometry of the channel section.
  // 
  // In quasilinear form: M * du/dt + (dQ/dA * dA/du) * du/dx = 0
  // C:= dQ/dA: the jacobian dQ/dA is the phase velocity or advective
  //                  jacobian (a scalar)
  // a:= the jacobian `dQ/du = dQ/dA * dA/du = C * M' is the
  //                  jacobian
  // M= `dA/du = wl_width:= ' is the analogous to `mass' or specific heat.
  double Q,CC,a,S;

  options |= SCALAR_TAU;	// tell the advective element routine
				// that we are returning a scalar tau
  // calls the `channel_shape' object
  set_state(U,grad_U);
  // slope is computed from the gradient of the bottom height (entered
  // as third component in the `nodes' section
  S = -grad_H.get(1,1);
  assert(S>=0.); // slope should be negative in the direction
		 // of the stream
  // computes `Q' and `C'
  friction_law->flow(area,perimeter,S,Q,C);
  a = C * wl_width;
  flux.ir(1,1).set(Q).rs();
  fluxd.set(0.);
  grad_U.ir(2,1);
  A_grad_U.set(grad_U).scale(a);
  grad_U.rs();

  if ( options & COMP_UPWIND ) {
    double h_supg = 2./iJaco.get(1,1);
    CC = fabs(C);
    tau_supg.setel(h_supg/(2.*CC),1,1);
    lam_max = CC;
  }

  if (options & COMP_SOURCE)  G_source.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0); // fixme:= Not implemented yet, used for
	     // absorbing boundary conditions
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  grad_N.ir(1,1);
  A_grad_N.ir(2,1).ir(3,1).set(grad_N).scale(C*wl_width);
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
ChannelShape *ChannelShape::factory(const NewElemset *elemset) {
  int ierr;
  //o Choose channel section shape, may be {\tt circular}
  // or {\tt rectangular}
  EGETOPTDEF(elemset,string,shape,string("undefined"));
  if (shape=="rectangular" || shape=="undefined")  {
    return new rect_channel(elemset);
  } else if (shape=="circular") {
    return new circular_channel(elemset);
  } else assert(0); // Not valid shape
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FrictionLaw::FrictionLaw(const NewElemset *e) : elemset(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FrictionLaw *FrictionLaw::factory(const NewElemset *elemset) {
  int ierr;
  //o Choose friction law, may be {\tt manning} or {\tt chezy}
  EGETOPTDEF(elemset,string,friction_law,string("undefined"));
  if (friction_law=="Chezy" || friction_law=="undefined")  {
    return new Chezy(elemset);
  } else if (friction_law=="Manning") {
    return new Manning(elemset);
  } else assert(0); // Not valid friction law
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Manning::init() { 
  int ierr;
  //o Roughness coefficient for the Manning formula (a.k.a. $n$)
  GETOPT_PROP(double,roughness,1.);
  // elemset->get_prop(roughness_prop,"roughness");
  //o Unit conversion factor for Manning friction law.
  EGETOPTDEF_ND(elemset,double,a_bar,1.);
  assert(ierr==0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void StreamLossFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
			  FastMat2 &jacin,FastMat2 &jacout) {
  flux.set(uout).rest(H_in).rest(uin).scale(k);
  jacin.set(k);
  jacout.set(-k);
}

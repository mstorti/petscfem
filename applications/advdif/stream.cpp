//__INSERT_LICENSE__

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "stream.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void DummyEnthalpyFun
::set_state(const FastMat2 &U) { s->UU.set(U); }

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
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::start_chunk(int &ret_options) {
  if (!channel) channel = ChannelShape::factory(elemset);
  channel->init();
  elemset->get_prop(slope_prop,"slope");

#if 0
  /// Friction related stuff
  string f = "Chezy";
  elemset->get_string("friction_law",f,1);
  if (f=="Chezy") {
    friction_law = Chezy;
    elemset->get_prop(Ch_prop,"Ch");
  } else if (f=="Manning") {
    friction_law = Manning;
    elemset->get_prop(roughness_prop,"roughness");
  } else assert(0);
#endif
  friction_law->init();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::element_hook(ElementIterator &element) {
  S = elemset->prop_val(element,slope_prop);
#if 0
  if (friction_law==Chezy) {
    Ch = elemset->prop_val(element,Ch_prop);
  } else {
    Ch = elemset->prop_val(element,Ch_prop);
  }
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stream_ff::stream_ff(const NewAdvDif *e) 
  : AdvDifFFWEnth(e), channel(NULL), 
    friction_law(FrictionLaw::factory(elemset)) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stream_ff::~stream_ff() { delete friction_law; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stream_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  u = U.get(1);
  channel->geometry(u,area,wl_width,perimeter);
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

  // Q:= the volumetric flow, and C:= the advective jacobian
  double Q,C;

  options |= SCALAR_TAU;	// tell the advective element routine
				// that we are returning a scalar tau
  // calls the `channel_shape' object
  set_state(U,grad_U);
  // computes `Q' and `C'
  friction_law->flow(area,perimeter,S,Q,C);
  flux.ir(1,1).set(Q).rs();
  fluxd.set(0.);
  grad_U.ir(2,1);
  A_grad_U.set(grad_U).scale(C);
  grad_U.rs();

  if ( options & COMP_UPWIND ) {
    double h_supg = 2./iJaco.get(1,1);
    tau_supg.setel(h_supg/(2.*C),1,1);
    lam_max = C;
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
  A_grad_N.ir(2,1).ir(3,1).set(grad_N).scale(wl_width);
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
ChannelShape *ChannelShape::factory(const NewElemset *e) {
  return new rect_channel(e);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void rect_channel::init() {
  elemset->get_prop(width_prop,"width");
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void rect_channel
::element_hook(const NewElemset *elemset,ElementIterator element) {
  width = elemset->prop_val(element,width_prop);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void rect_channel
::geometry(double u,double &area,double &wl_width,double &perimeter) {
  area = u*width;
  wl_width = width;
  perimeter = width+2*u;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FrictionLaw::FrictionLaw(const NewElemset *e) : elemset(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
FrictionLaw *FrictionLaw::factory(const NewElemset *e) {
  return new Chezy(e);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Chezy::init() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Chezy::element_hook(const NewElemset *elemset,
			 ElementIterator element) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Chezy::flow(double area,double perimeter,double S,
			double &Q,double &C) const {
  double Ch = 1., m=1.5;
  double gamma = Ch*sqrt(S/perimeter);
  Q = gamma*pow(area,m);
  C = Q*m/area;
}

//__INSERT_LICENSE__
//$Id: svenant1d.cpp,v 1.1 2005/08/24 01:52:49 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/generror.h>

#include "svenant1d.h"

extern int MY_RANK,SIZE;

#define GF_GETOPTDEF_ND(type,var,def)				\
 { if (elemset) { EGETOPTDEF_ND(elemset,type,var,def); }	\
 else if (old_elemset)						\
      { TGETOPTDEF_ND(old_elemset->thash,type,var,def); }	\
 else assert(0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::start_chunk(int &ret_options) {
  int ierr;

  if (elemset) {
    new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
    elemset->elem_params(nel,ndof,nelprops);
  } else if (old_elemset) {
    nel = old_elemset->nel;
    ndof = old_elemset->ndof;
  } else assert(0);

  GF_GETOPTDEF_ND(int,ndim,0);

  double channel_width;
  //o Width of channel
  GF_GETOPTDEF_ND(double,channel_width,0.0);
  assert(channel_width>0.);
  w = channel_width;

  //o Gravity acceleration
  GF_GETOPTDEF_ND(double,gravity,0.0);
  assert(gravity>0.);

  assert(ndim==1);
  assert(ndof==2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
svenant1d_ff::svenant1d_ff(NewElemset *e)
  : AdvDifFFWEnth(e), old_elemset(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
svenant1d_ff::svenant1d_ff(Elemset *e)
  : old_elemset(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
svenant1d_ff::~svenant1d_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::set_state(const FastMat2 &UU) {
  U.set(UU);
  double dummy;
  // Scalar variables
  A = U.get(1);
  Q = U.get(2);
  u = Q/A;

  // This is specific for rectangular channel
  h = A/w;
  F = w*h*h/2;
  wl = w;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
set_state(const FastMat2 &U,
	  const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::enthalpy(FastMat2 &H) {
  H.set(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,
	      const FastMat2 &N,double weight) {
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
comp_P_Cp(FastMat2 &P_Cp,
	  const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
get_Cp(FastMat2 &Cp_a) {
  Cp_a.set(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
get_Ajac(FastMat2 &Ajac_a) {
  Ajac_a.set(Ajac);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#if 0
void svenant1d_ff::
compute_tau(int ijob,double &delta_sc) {
  const FastMat2 &grad_N = *advdf_e->grad_N();
  double tol=1.0e-16;
  velmod = Q/A;
  h_supg = grad_N.sum_abs_all();
  delta_sc = 0.5*h_supg*velmod;
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  Cp.eye();
  flux.setel(Q,1);
  flux.setel(Q*Q/A+gravity*F,2);

  Ajac.setel(0.0,1,1,1);
  Ajac.setel(1.0,1,1,2);
  Ajac.setel(u*Q+gravity*A/wl,1,2,1);
  Ajac.setel(2*u,1,2,2);

  fluxd.set(0.0);

  if (options & COMP_SOURCE) {
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(Ajac,normal,-1,1,2,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(Ajac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				     FastMat2 &dshapex,double w) {
  grad_N_D_grad_N.set(0.0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  N_N_C.set(0.0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  N_P_C.set(0.0);
}

#undef SHV
#define SHV(mess,v) printf("[%d] %s " #v " %f\n",	\
	    MY_RANK,mess,v.sum_square_all());

#ifdef USE_COMP_P_SUPG
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::comp_P_supg(FastMat2 &P_supg) {


}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
	    FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
compute_shock_cap_aniso(double &delta_aniso,
			FastMat2 &jvec_a) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void svenant1d_ff::
get_C(FastMat2 &C) {
  C.set(0.);
}

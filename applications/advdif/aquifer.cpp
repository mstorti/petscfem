//__INSERT_LICENSE__
//$Id: aquifer.cpp,v 1.18 2005/05/02 00:36:49 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/generror.h>

#include "aquifer.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::start_chunk() {
  int ierr;
  
  EGETOPTDEF_ND(elemset,int,ndim,0); //nd
  //o Constant rain
  EGETOPTDEF_ND(elemset,double,rain,0.); 
  //o Threshold for wet aquifer width ( #phi-eta# )
  EGETOPTDEF_ND(elemset,double,wet_aquifer_width_min,0.); 
  assert(wet_aquifer_width_min>=0.);
  //o Flag ehether to stop on dry aquifer
  EGETOPTDEF_ND(elemset,int,dry_aquifer_stop,0); 
  //o Threshold for wet aquifer width ( #phi-eta# )
  EGETOPTDEF_ND(elemset,double,wet_aquifer_width_min,0.); 
  assert(wet_aquifer_width_min>=0.);
  //o Flag ehether to stop on dry aquifer
  EGETOPTDEF_ND(elemset,int,dry_aquifer_stop,0); 
  //o We can have the #phi# dof of the aquifer to be
  //  any of the dofs there. This indx flags which of the dof's
  //  of the aquifer is the #phi# (piezometric height) field. 
  EGETOPTDEF_ND(elemset,int,aquifer_phi_dof,1); 
  assert(ndim==2);

  elemset->elem_params(nel,ndof,nelprops);
  assert(ndof>=aquifer_phi_dof);

  // elemset->get_prop(eta_pr,"eta");
  elemset->get_prop(K_pr,"K");
  elemset->get_prop(S_pr,"S");
  elemset->get_prop(rain_pr,"rain");

  tmp.resize(2,ndim,nel);
  tmp1.resize(2,nel,nel);

  tmp2.resize(1,nel);
  tmp3.resize(2,nel,nel);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::end_chunk() {
  tmp.clear();
  tmp1.clear();
  tmp2.clear();
  tmp3.clear();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
element_hook(ElementIterator &element) {
  S = elemset->prop_val(element,S_pr);
  K = elemset->prop_val(element,K_pr);
  element_m = element;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
gp_hook(int ipg,const FastMat2 &U,
	const FastMat2 &grad_U) {
  rain = 
    elemset->prop_val(element_m,
		      rain_pr,elemset->time());
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
compute_flux(const FastMat2 &U,const FastMat2
	     &grad_U,FastMat2 &fluxd,FastMat2 &G,
	     FastMat2 &H, FastMat2 &grad_H) {
  eta = H.get(1);
  phi = U.get(aquifer_phi_dof);
  ///  phi = U.get(1);
  double wet_aquifer_width = phi-eta;
  if(dry_aquifer_stop && phi-eta<=0.) 
    throw GenericError("Aquifers got dried.");
  if (wet_aquifer_width<wet_aquifer_width_min) 
    wet_aquifer_width = wet_aquifer_width_min;
  fluxd.set(grad_U).scale(K*wet_aquifer_width);
  G.set(rain);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		     FastMat2 & dshapex,double w) {
  tmp.set(dshapex).scale(w*K*(phi-eta));
  int d = aquifer_phi_dof;
  grad_N_D_grad_N.set(0.).ir(2,d).ir(4,d);
  grad_N_D_grad_N.prod(tmp,dshapex,-1,1,-1,2);
  grad_N_D_grad_N.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
enthalpy(FastMat2 &H, FastMat2 &U) {
  H.set(U).scale(S*(phi-eta));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::
comp_N_Cp_N(FastMat2 &N_Cp_N,
	    FastMat2 &N, double w) {
  tmp2.set(N).scale(S*(phi-eta)*w);
  tmp3.prod(N,tmp2,1,2);
  int d = aquifer_phi_dof;
  N_Cp_N.set(0.).ir(2,d).ir(4,d).set(tmp3);
  N_Cp_N.rs();
}

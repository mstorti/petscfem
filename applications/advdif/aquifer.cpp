//__INSERT_LICENSE__
//$Id: aquifer.cpp,v 1.7 2002/04/12 14:44:17 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "aquifer.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::start_chunk() {
  int ierr;
  
  EGETOPTDEF_ND(elemset,int,ndim,0); //nd
  assert(ndim==2);

  elemset->elem_params(nel,ndof,nelprops);
  assert(ndof==1);

  elemset->get_prop(eta_pr,"eta");
  elemset->get_prop(K_pr,  "K");
  elemset->get_prop(S_pr,  "S");

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
void aquifer_ff::element_hook(ElementIterator &element) {
  S = elemset->prop_val(element,S_pr);
  K = elemset->prop_val(element,K_pr);
  eta = elemset->prop_val(element,eta_pr);
#if 0
  int p,q;
  element.position(p,q);
  printf("element: %d, eta %f\n",p,eta);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::gp_hook(int ipg,const FastMat2 &U,const FastMat2 &grad_U) {
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::compute_flux(const FastMat2 &U,const FastMat2
			      &grad_U,FastMat2 &fluxd,FastMat2 &G,
			      FastMat2 &H, FastMat2 &grad_H) {
  eta = H.get(1);
  phi = U.get(1);
  assert(phi>=eta);
  fluxd.set(grad_U).scale(K*(phi-eta));
  G.set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp.set(dshapex).scale(w*K*(phi-eta));
  grad_N_D_grad_N.ir(2,1).ir(4,1);
  grad_N_D_grad_N.prod(tmp,dshapex,-1,1,-1,2);
  grad_N_D_grad_N.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::enthalpy(FastMat2 &H, FastMat2 &U) {
  H.set(U).scale(S*(phi-eta));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void aquifer_ff::comp_N_Cp_N(FastMat2 &N_Cp_N,FastMat2 &N, double w) {
  tmp2.set(N).scale(S*(phi-eta)*w);
  tmp3.prod(N,tmp2,1,2);
  N_Cp_N.ir(2,1).ir(4,1).set(tmp3);
  N_Cp_N.rs();
}


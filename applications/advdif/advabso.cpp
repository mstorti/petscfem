//__INSERT_LICENSE__
// $Id: advabso.cpp,v 1.1 2005/01/26 18:39:31 mstorti Exp $
#include "./advabso.h"

#define gasflow_abso gasflow_abso2

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
lm_initialize() { 
  int ff_options=0;
  adv_diff_ff->start_chunk(ff_options);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
init() {
  int ierr;
  TGETOPTDEF_ND(thash,int,ndim,0);
  flux.resize(2,ndof,ndim);
  fluxd.resize(2,ndof,ndim);
  A_grad_U.resize(1,ndof);
  grad_U.resize(2,ndim,ndof).set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
res(int k,FastMat2 &U,FastMat2 &r,
    FastMat2 &w,FastMat2 &jac) {
  double delta_sc=0.0,
      lambda_max_pg=0.0;
  U.ir(1,1);
  Ucpy.set(U);
  U.rs();
  adv_diff_ff->set_state(Ucpy,grad_U);
  adv_diff_ff
    ->compute_flux(Ucpy, dummy, dummy, dummy, flux, fluxd,
		   A_grad_U, grad_U, dummy,
		   dummy, delta_sc, lambda_max_pg, dummy,
		   dummy, dummy, dummy, 0);
}

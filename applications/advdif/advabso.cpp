//__INSERT_LICENSE__
// $Id: advabso.cpp,v 1.3 2005/01/27 14:11:34 mstorti Exp $
#include "./advabso.h"

#define gasflow_abso gasflow_abso2

static
double p(double x) {
  return (x>=0.? x : 0.);
}

static
double m(double x) {
  return (x<0. ? x : 0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
lm_initialize() { 
  int ff_options=0;
  adv_diff_ff->start_chunk(ff_options);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
init() {
  int nelprops;
  elem_params(nel,ndof,nelprops);
  assert(nel==3);
  int ierr;
  NSGETOPTDEF_ND(int,ndim,0);
  flux.resize(2,ndof,ndim);
  fluxd.resize(2,ndof,ndim);
  A_grad_U.resize(1,ndof);
  grad_U.resize(2,ndim,ndof).set(0.);
  normal.resize(1,ndim);
  A_jac.resize(2,ndof,ndof);
  S.resize(2,ndof,ndof);
  invS.resize(2,ndof,ndof);
  tmp1.resize(2,ndof,ndof);
  c.resize(2,2,ndof);
  cm.resize(1,ndof);
  cp.resize(1,ndof);
  Pi_m.resize(2,ndof,ndof);
  Pi_p.resize(2,ndof,ndof);
  get_prop(normal_prop,"normal");
  // The state on the reference node
  Uref.resize(1,ndof);
  dU.resize(1,ndof);
  // The lagrange multipliers (state for the 2nd node)
  Ulambda.resize(1,ndof);
  // The state of the outlet node
  Uo.resize(1,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
res(int k,FastMat2 &U,FastMat2 &r,
    FastMat2 &w,FastMat2 &jac) {
  double delta_sc=0.0,
      lambda_max_pg=0.0;
  U.ir(1,1); Uo.set(U);
  U.ir(1,2); Ulambda.set(U);
  U.ir(1,3); Uref.set(U);
  U.rs();
  adv_diff_ff->set_state(Uref,grad_U);
  adv_diff_ff
    ->compute_flux(Ucpy, dummy, dummy, dummy, flux, fluxd,
		   A_grad_U, grad_U, dummy,
		   dummy, delta_sc, lambda_max_pg, dummy,
		   dummy, dummy, dummy, 0);
  adv_diff_ff->comp_A_jac_n(A_jac,normal);
  c.eig(A_jac,S);
  invS.inv(S);
  double aimag = c.ir(1,2).sum_square_all();
  assert(aimag==0.0);
  c.ir(1,1);
  // Pi_m = projector on negative eigenvalues space
  // Pi_p = projector on positive eigenvalues space
  // Pi_m + Pi_p = A_jac
  Pi_m.set(0.).d(1,2).set(c).rs();
  Pi_p.set(Pi_m);
  Pi_m.fun(m);
  Pi_p.fun(p);
  tmp1.prod(Pi_m,invS,1,-1,-1,2);
  Pi_m.prod(S,tmp1,1,-1,-1,2);
  tmp1.prod(Pi_p,invS,1,-1,-1,2);
  Pi_p.prod(S,tmp1,1,-1,-1,2);
  tmp1.set(Pi_m).add(Pi_p);
  // residual is the projection of U-Uref
  // on to the space of incoming waves
  dU.set(Uo).rest(Uref);
  r.prod(Pi_m,dU,1,-1,-1);
  // The vector of reactions is the pojector on
  // to the incoming wave space
  w.set(0.).ir(1,1).set(Pi_m).rs();
  jac.ir(2,1).set(Pi_m).rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
element_hook(ElementIterator &element) {
  assert(normal_prop.length == ndim);
  normal.set(prop_array(element,normal_prop));
}

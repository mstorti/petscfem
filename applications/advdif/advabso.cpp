//__INSERT_LICENSE__
// $Id: advabso.cpp,v 1.8 2005/01/28 18:16:44 mstorti Exp $
#include "./advabso.h"

#define gasflow_abso gasflow_abso2

// static
// double p(double x) {
//   return (x>=0.? x : 0.);
// }

static
double msign(double x) {
  return (x<0. ?  1. : 0.0);
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
  //o Flags whether to use the old state ad the
  // boundary as reference for the linear absorbing
  // boundary condition.
  NSGETOPTDEF_ND(int,use_old_state_as_ref,0);
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
  Pi_m.resize(2,ndof,ndof);
  Pi_p.resize(2,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  Uold.resize(2,nel,ndof);
  invCp.resize(2,ndof,ndof);
  get_prop(normal_prop,"normal");
  assert(normal_prop.length == ndim);
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
  if (use_old_state_as_ref) {
    get_old_state(Uold);
    Uold.ir(1,1);
    Uref.set(Uold);
    Uold.rs();
  } else {
    U.ir(1,3); Uref.set(U);
  }
  U.rs();
  adv_diff_ff->set_state(Uref,grad_U);
  adv_diff_ff
    ->compute_flux(Uref, dummy, dummy, dummy, flux, fluxd,
		   A_grad_U, grad_U, dummy,
		   dummy, delta_sc, lambda_max_pg, dummy,
		   dummy, dummy, dummy, 0);
  adv_diff_ff->comp_A_jac_n(A_jac,normal);
  adv_diff_ff->get_Cp(Cp);
  invCp.inv(Cp);
  // tmp1 = Cp \ A
  tmp1.prod(invCp,A_jac,1,-1,-1,2);
  c.eig(tmp1,S);
  invS.inv(S);
  double aimag = c.ir(1,2).sum_square_all();
  assert(aimag<1e-10);
  c.ir(1,1);
  // Pi_m = projector on negative eigenvalues space
  // Pi_p = projector on positive eigenvalues space
  // Pi_m + Pi_p = A_jac
  Pi_m.set(0.).d(1,2)
    .set(c).fun(msign).rs();
  c.rs();
  tmp1.prod(Pi_m,invS,1,-1,-1,2);
  Pi_m.prod(S,tmp1,1,-1,-1,2);
#if 0
  tmp1.prod(Pi_p,invS,1,-1,-1,2);
  Pi_p.prod(S,tmp1,1,-1,-1,2);
  tmp1.set(Pi_m).add(Pi_p);
#endif
  // residual is the projection of U-Uref
  // on to the space of incoming waves
  dU.set(Uo).rest(Uref);
  r.prod(Pi_m,dU,1,-1,-1);
  // The vector of reactions is the pojector on
  // to the incoming wave space: w = Cp * Pi_m
  tmp1.prod(Cp,Pi_m,1,-1,-1,2);
  w.set(0.).ir(1,1).set(tmp1).rs();
  jac.ir(2,1).set(Pi_m);
  if (!use_old_state_as_ref) 
    jac.ir(2,3).set(Pi_m).scale(-1.0);
  jac.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void AdvectiveAbso::
element_hook(ElementIterator &element) {
  normal.set(prop_array(element,normal_prop));
}

//__INSERT_LICENSE__
//$Id: ffburg.cpp,v 1.16 2002/01/14 03:45:05 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "burgers.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void UGlobalAJac::init(int ndim_,int ndof_,int nel_) {
  ndim = ndim_;
  ndof = ndof_;
  nel = nel_;
  tmp2.resize(1,ndof).set(1.);
}

void UGlobalAJac::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(U,u,1,2);
}

void UGlobalAJac::
comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  tmp5.prod(u,normal,-1,-1);
  double un = tmp5.get();
  A_jac_n.eye(un);
}

void UGlobalAJac::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(u,grad_U,-1,-1,1);
}
  
void UGlobalAJac::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  tmp.prod(u,dshapex,-1,-1,1);
  for (int j=1; j<=nel; j++) {
    A_grad_N.ir(1,j).eye(tmp.get(j));
  }
  A_grad_N.rs();
}
  
void UGlobalAJac::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  tmp3.prod(iJaco,u,1,-1,-1);
  Uintri.prod(tmp2,tmp3,1,2);
}

void UGlobalAJac::
comp_vel_per_field(FastMat2 &vel_per_field) {
  double vv= u.sum_square_all();
  vel_per_field.set(sqrt(vv));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void GlobalScalarDJac::init(int ndim_,int ndof_,int nel_) {
  ndim = ndim_;
  ndof = ndof_;
  nel = nel_;
  eye_ndof.resize(2,ndof,ndof).eye(1.);
}

void GlobalScalarDJac
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).scale(diff).rs();
}

void GlobalScalarDJac
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp.prod(dshapex,dshapex,-1,1,-1,2).scale(w*diff);
  grad_N_D_grad_N.prod(tmp,eye_ndof,1,3,2,4);
}

void GlobalScalarDJac
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(diff);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
advdif_wjac_ff::advdif_wjac_ff(NewElemset *elemset_,
			       DJac *d=NULL,AJac *a=NULL, 
			       CJac *c=NULL, SourceTerm *st=NULL) 
  : NewAdvDifFF(elemset_), a_jac(a),
    source_term(st), c_jac(c), d_jac(d) {};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#if 0
#undef __FUNC__
#define __FUNC__ "burgers_ff::burgers_ff()"
burgers_ff::burgers_ff(NewElemset *elemset_) 
  : advdif_wjac_ff(elemset_), a_jac(&u_a_jac),
    source_term(&st), c_jac(&n_c_jac), d_jac(&g_d_jac) {};
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void burgers_ff::start_chunk(int &ret_options) {
  
  FastMat2 &u = u_a_jac.u;
  int ierr;
  EGETOPTDEF_ND(elemset,int,ndim,0); //nd

  //o Scale the SUPG upwind term. Set to 0 in order to
  //  not to include the upwind term. 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  elemset->elem_params(nel,ndof,nelprops);

  enthalpy_fun = &identity_ef;
  identity_ef.init(ndof,ndim,nel);

  u_a_jac.init(ndim,ndof,nel);
  g_d_jac.init(ndim,ndof,nel);
  // this must be done here. Because this matrices may have
  // differente sizes for different elemsets. 
  u0.resize(1,ndim);
  u.resize(1,ndim);
  Uintri.resize(1,ndim);
  tmp0.resize(1,ndim);

  assert(ndof==1); // Only 1D Burgers considered

  //o Diffusivity (viscosity)
  EGETOPTDEF_ND(elemset,double,diffusivity,0.);
  g_d_jac.diff = diffusivity;

  //o Flux law is $f= 0.5\,c\, \phi^2\,!u_0$, where $c$
  // is this coefficient. 
  EGETOPTDEF_ND(elemset,double,flux_law_coefficient,1.);
  //o Scale the SUPG upwind term. 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  //o _T: double[ndim] _N: u0 _D: [1,0] 
  // _DOC: Vector defining direction for flux. 
  // _END
  u0.set(0.); u0.setel(1.,1);
  ierr = elemset->get_double("u0",*u0.storage_begin(),1,ndim); 
  PETSCFEM_ASSERT0(ierr==0,"Error reading variable u0");
  ret_options |= SCALAR_TAU;	// tell the advective element routine
				// that we are returning a scalar tau
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void burgers_ff::compute_flux(...)"
void burgers_ff::compute_flux(COMPUTE_FLUX_ARGS) {
  FastMat2 &u = u_a_jac.u;

  double coef = flux_law_coefficient;
  double phi = U.get(1);
  u.set(u0).scale(coef*phi);
  double vel = sqrt(u.sum_square_all());
  Ucpy.set(U);

  // A_jac.ir(2,1).ir(3,1).set(u).rs();
  flux.ir(1,1).set(u0).scale(0.5*coef*phi*phi).rs();

  if (options & COMP_UPWIND) {

#if 0
    D_jac.set(0.).ir(3,1).ir(4,1).eye(diffusivity).rs();
    grad_U.t();
    fluxd.set(grad_U).scale(diffusivity);
    grad_U.rs();
    A_grad_U.prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);
#endif

    d_jac->comp_fluxd(fluxd,grad_U);
    a_jac->comp_A_grad_U(A_grad_U,grad_U);

    lam_max = vel;
    tau_supg.set(0.);

    Uintri.prod(iJaco,u,1,-1,-1);
    double Uh = sqrt(Uintri.sum_square_all()); // this is
				// approx. 2*U/h

    double &alpha = diffusivity;
    double tau;

    FastMat2::branch();
    if (vel*vel > 1e-5*Uh*alpha) { // remove singularity when v=0
      FastMat2::choose(0);
      double Pe  = vel*vel/(Uh*alpha);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = tau_fac/Uh*magic; // intrinsic time
    } else {
      FastMat2::choose(1);
      double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
      tau = tau_fac*h*h/(12.*alpha);
    }
    FastMat2::leave();
    tau_supg.setel(tau,1,1);
    delta_sc = 0.;
  }
  
  if (options & COMP_SOURCE) {
    G_source.set(0.);		// Only null source term allowed
				// right now!!
  }
}

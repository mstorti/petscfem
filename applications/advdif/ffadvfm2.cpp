/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include "../../src/fem.h"
#include "../../src/texthash.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"
#include "../../src/util2.h"
#include "../../src/fastmat2.h"

#include "nwadvdif.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
newadvecfm2_ff_t::newadvecfm2_ff_t(NewAdvDif *elemset_) 
  : NewAdvDifFF(elemset_), u_per_field(*this), u_global(*this), 
  full_adv_jac(*this), full_dif_jac(*this),
  scalar_dif_per_field(*this), global_scalar_djac(*this)
{};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::GlobalScalar
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.t().set(grad_U).scale(*(ff.difjac)).rs();
}

void newadvecfm2_ff_t::GlobalScalar
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  tmp.prod(dshapex,dshapex,-1,1,-1,2).scale(w*(*(ff.difjac)));
  grad_N_D_grad_N.prod(tmp,ff.eye_ndof,1,3,2,4);
#if 0
  ff.tmp3.set(*(ff.difjac));
  D_grad_N.d(3,2).prod(dshapex,ff.tmp3,1,3,2).rs();
#endif
}

void newadvecfm2_ff_t::GlobalScalar
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(*(ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::ScalarDifPerField
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.set(grad_U);
  for (int j=0; j<ff.elemset->ndof; j++) 
    fluxd.ir(1,j).scale(ff.D_jac.get(j));
  fluxd.rs();
}

void newadvecfm2_ff_t::ScalarDifPerField
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
#if 0
  D_grad_N.prod(ff.D_jac,dshapex,1,-1,2,3,-1,4);
#endif
}

void newadvecfm2_ff_t::ScalarDifPerField
::comp_dif_per_field(FastMat2 &dif_per_field) {
  dif_per_field.set(*(ff.difjac));
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::FullDifJac
::comp_fluxd(FastMat2 &fluxd,FastMat2 &grad_U) {
  fluxd.prod(ff.D_jac,grad_U,2,-1,1,-2,-1,-2);
}

void newadvecfm2_ff_t::FullDifJac
::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
		       FastMat2 & dshapex,double w) {
  D_grad_N.prod(ff.D_jac,dshapex,1,-1,2,3,-1,4)
    .scale(w);
  grad_N_D_grad_N.prod(D_grad_N,dshapex,-1,2,4,1,-1,3);
}

void newadvecfm2_ff_t::FullDifJac
::comp_dif_per_field(FastMat2 &dif_per_field) {
  ff.D_jac.d(4,3);
  dif_per_field.sum_abs(ff.D_jac,-1,-2,1).scale(1./ff.elemset->ndim);
  ff.D_jac.rs();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::FullAdvJac::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(ff.u,U,2,1,-1,-1);
}

void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,1,-2,-1,-2);
}
  
void newadvecfm2_ff_t::FullAdvJac::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.prod(dshapex,ff.u,-1,1,-1,2,3);
}

void newadvecfm2_ff_t::FullAdvJac::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  // Here we take the diagonal of each Jacobian as the component
  // of a per field vector velocity
  ff.u.d(3,2);
  Uintri.prod(iJaco,ff.u,2,-1,-1,1);
  ff.u.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::UPerField::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.set(ff.u);
  for (int j=1; j<=ff.elemset->ndof; j++) {
    flux.ir(1,j).scale(U.get(j));
  }
  flux.rs();
}

void newadvecfm2_ff_t::UPerField::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  for (int j=1; j<=ff.elemset->ndof; j++) {
    grad_U.ir(2,j);
    ff.u.ir(1,j);
    tmp.prod(grad_U,ff.u,-1,-1);
    A_grad_U.setel(tmp.get(),j);
  }
  grad_U.rs();
  ff.u.rs();
  A_grad_U.rs();
}
  
void newadvecfm2_ff_t::UPerField::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  A_grad_N.set(0.).d(3,2).prod(dshapex,ff.u,-1,1,2,-1).rs();
}

void newadvecfm2_ff_t::UPerField::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  Uintri.prod(iJaco,ff.u,2,-1,1,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:   
void newadvecfm2_ff_t::UGlobal::comp_flux(FastMat2 &flux,FastMat2 &U) {
  flux.prod(U,ff.u,1,2);
}

void newadvecfm2_ff_t::UGlobal::
comp_A_grad_U(FastMat2 &A_grad_U,FastMat2 &grad_U) {
  A_grad_U.prod(ff.u,grad_U,-1,-1,1);
}
  
void newadvecfm2_ff_t::UGlobal::
comp_A_grad_N(FastMat2 &A_grad_N,FastMat2 &dshapex) {
  tmp.prod(ff.u,dshapex,-1,-1,1);
  for (int j=1; j<=ff.elemset->nel; j++) {
    A_grad_N.ir(1,j).eye(tmp.get(j));
  }
  A_grad_N.rs();
}
  
void newadvecfm2_ff_t::UGlobal::
comp_Uintri(FastMat2 &Uintri,FastMat2 &iJaco) {
  tmp3.prod(iJaco,ff.u,1,-1,-1);
  Uintri.prod(ff.tmp2,tmp3,1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::element_hook(ElementIterator &)"
// This is to pass to the advective function the element
void newadvecfm2_ff_t::element_hook(ElementIterator &element_) {
  element = element_;
  advjac = elemset->prop_array(element,advective_jacobians_prop);
  u.set(advjac);

  difjac = elemset->prop_array(element,diffusive_jacobians_prop);
  d_jac->update(difjac);
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void advecfm2_ff_t::start_chunk(int ret_options)"
void newadvecfm2_ff_t::start_chunk(int ret_options) {
  // FastMat2Shell *A_jac_l = elemset->A_jac;
  ndim = elemset->ndim;
  ndof = elemset->ndof;
  Uintri.resize(2,ndof,ndim);
  tmp2.resize(1,ndof).set(1.);
  tmp3.set(tmp2);
  dif_per_field.resize(1,ndof);

  ret_options &= !SCALAR_TAU; // tell the advective element routine

  // Read advective jacobians
  //o _T: double[ndim]/double[ndim*ndof]/double[ndim*ndof*ndof] 
  //  _N: advective_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex advective_jacobians
  //  _END
  elemset->get_prop(advective_jacobians_prop,"advective_jacobians");

  if (advective_jacobians_prop.length == ndim) {
    u.resize(1,ndim);
    a_jac =  &u_global;
  } else if (advective_jacobians_prop.length == ndim*ndof) {
    u.resize(2,ndof,ndim);
    a_jac =  &u_per_field;
  } else if (advective_jacobians_prop.length == ndim*ndof*ndof) {
    u.resize(3,ndim,ndof,ndof);
    a_jac =  &full_adv_jac;
  } else {
    assert(0);
  }

  // Read diffusive jacobians (diffusivity matrices)
  //o _T: double[var_len] 
  //  _N: diffusive_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex diffusive_jacobians
  //  _END
  elemset->get_prop(diffusive_jacobians_prop,"diffusive_jacobians");

  if (diffusive_jacobians_prop.length == 1) {
    d_jac =  &global_scalar_djac;
    eye_ndof.resize(2,ndof,ndof).eye(1.);
  } else if (diffusive_jacobians_prop.length == ndof) {
    D_jac.resize(1,ndof);
    d_jac =  &scalar_dif_per_field;
  } else if (diffusive_jacobians_prop.length == ndim*ndim*ndof*ndof) {
    D_jac.resize(4,ndim,ndim,ndof,ndof);
    d_jac =  &full_dif_jac;
  } else {
    assert(0);
  }

  // Read reactive jacobians (reactive  matrix)
  //o _T: double[var_len] 
  //  _N: reactive_jacobians _D: all zero  _DOC: 
  //  FIXME:= TO BE DOCUMENTED LATER
  //  _END
  C_jac_l.resize(2,ndof,ndof);
  const char *reaje;
  VOID_IT(cjacv);
  elemset->get_entry("reactive_jacobians",reaje); 
  if (reaje==0) {
    for (int k=0; k<ndof; k++) cjacv.push_back(0.);
  } else {
    read_double_array(cjacv,reaje); 
  }
  cjacvp=cjacv.begin();
  C_jac_l.set(0.);

  //o Scale the SUPG upwind term. 
  int ierr;
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  nc = cjacv.size();

  if (nc==ndof) {
    ret_options &= !SCALAR_TAU; // tell the advective element routine
				// that we are returning a non-scalar tau
    C_jac_l.set(0.);
    for (int k=1; k<=ndof; k++) {
      double beta=cjacv[k-1];
      C_jac_l.setel(beta,k,k);
    }

  } else {
    assert(0); // Not implemented yet
  }

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void newadvecfm2_ff_t::compute_flux(COMPUTE_FLUX_ARGS) {
  int ierr;
  double tau_a, tau_delta, gU, A01v[9];

  // Unfortunately we have to use copies of U and
  // iJaco due to const'ness restrictions.

  Ucpy.set(U);
  iJaco_cpy.set(iJaco);
  a_jac->comp_flux(flux,Ucpy);

  if (options & COMP_UPWIND) {

    d_jac->comp_fluxd(fluxd,grad_U);

    // A_grad_U es ndof x 1
    // A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);
    a_jac->comp_A_grad_U(A_grad_U,grad_U);

    lam_max = 0.;
    tau_supg.set(0.);

    a_jac->comp_Uintri(Uintri,iJaco_cpy);
    d_jac->comp_dif_per_field(dif_per_field);

    for (int k=1; k<=ndof; k++) {
      
      Uintri.ir(1,k);
      double alpha=dif_per_field.get(k);
      double vel = sqrt(u.sum_square_all());

      // double h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
      double Uh = sqrt(Uintri.sum_square_all()); // this is
				// approx. 2*U/h
      double tau;
      FastMat2::branch();
      if (vel*vel > 1e-5*Uh*alpha) {		// remove singularity when v=0
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
      tau_supg.setel(tau,k,k); // Setting in the intrinsic time
				// matrix
      if (vel>lam_max) lam_max = vel;
    }
    Uintri.rs();
    delta_sc = 0.;
  }
  
  if (options & COMP_SOURCE) {
    G_source.set(0.);		// Only null source term allowed
				// right now!!

    C_jac.set(C_jac_l);
    G_source.prod(C_jac,U,1,-1,-1).scale(-1.);

  }
}

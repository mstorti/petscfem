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

#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "flux_fun_burgers" 
int flux_fun_burgers(AD_FLUX_FUN_ARGS) {

  int ierr;

  static double diffusivity, tau_fac, flux_law_coefficient;
  static int ndof;
  static FastMat2 u0,u,Uintri,tmp0;

  // Load properties only once.
  if (start_chunk) {
    FastMat2::deactivate_cache();
    start_chunk = 0;

    // this must be done here. Because this matrices may have
    // differente sizes for different elemsets. 
    u0.resize(1,ndim);
    u.resize(1,ndim);
    Uintri.resize(1,ndim);
    tmp0.resize(1,ndim);

    ndof = U.dim(1);
    assert(ndof==1); // Only 1D Burgers considered
    // assert(ndim==1);

    //o Diffusivity (viscosity)
    SGETOPTDEF_ND(double,diffusivity,0.);
    //o Flux law is $f= 0.5\,c\, \phi^2\,!u_0$, where $c$
    // is this coefficient. 
    SGETOPTDEF_ND(double,flux_law_coefficient,1.);
    //o Scale the SUPG upwind term. 
    SGETOPTDEF_ND(double,tau_fac,1.);

    //o _T: double[ndim] _N: u0 _D: unit vector along $x$ axis 
    // _DOC: Vector defining direction for flux. 
    // _END
    u0.set(0.); u0.setel(1.,1);
    ierr = get_double(thash,"u0",u0.storage_begin(),1,ndim); CHKERRQ(ierr);
    FastMat2::activate_cache();
  }

  ret_options |= SCALAR_TAU;	// tell the advective element routine
				// that we are returning a scalar tau
  double coef = flux_law_coefficient;
  double phi = U.get(1);
  u.set(u0).scale(coef*phi);
  double vel = sqrt(u.sum_square_all());

  A_jac.ir(2,1).ir(3,1).set(u).rs();
  flux.ir(1,1).set(u0).scale(0.5*coef*phi*phi).rs();

  if (options & COMP_UPWIND) {

    D_jac.set(0.).ir(3,1).ir(4,1).eye(diffusivity).rs();
    grad_U.t();
    fluxd.set(grad_U).scale(diffusivity);
    grad_U.rs();

    A_grad_U.prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);

    lam_max = vel;
    tau_supg.set(0.);

    Uintri.prod(iJaco,u,1,-1,-1);
    double Uh = sqrt(Uintri.sum_square_all()); // this is
				// approx. 2*U/h

    double alpha = diffusivity;
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

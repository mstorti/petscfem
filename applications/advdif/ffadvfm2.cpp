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
#define __FUNC__ "advecfm2_ff_t::operator()" 
int advecfm2_ff_t::operator()(ADVDIFFF_ARGS) {

  static double ajacx[9],ajacy[9];
  int ierr;

  static int flag=0,ndof;
  double tau_a, tau_delta, gU, A01v[9];
  static int shock_capturing,na,nd;
  static FastMat2 A_jac_l, D_jac_l, tmp0;
  static double tau_fac;
  static FastMat2 u,u2,Uintri(1,ndim),AA;
  static vector<double> ajacv,djacv;
  static double *ajacvp,*djacvp;

  // Load properties only once.
  FastMat2::branch();
  if (start_chunk) {
    FastMat2::choose(0);
    start_chunk = 0;

    ndof = U.dim(1);
    // Read advective jacobians
    A_jac_l.resize(3,ndim,ndof,ndof);
    u.resize(2,ndof,ndim);
    //o _T: double[ndim]/double[ndim*ndof]/double[ndim*ndof*ndof] 
    //  _N: advective_jacobians _D: no default  _DOC: 
    //i_tex advdif.tex advective_jacobians
    //  _END
    const char *advje;
    VOID_IT(ajacv);
    elemset->get_entry("advective_jacobians",advje); CHKERRQ(advje==0);
    read_double_array(ajacv,advje);
    ajacvp=ajacv.begin();

    // Read diffusive jacobians (diffusivity matrices)
    D_jac_l.resize(4,ndim,ndim,ndof,ndof);
    //o _T: double[ndof]/double[ndim*ndim*ndof]/double[ndim*ndim*ndof*ndof] 
    //  _N: diffusive_jacobians _D: no default  _DOC: 
    //i_tex advdif.tex advective_jacobians
    //  _END
    const char *difje;
    VOID_IT(djacv);
    elemset->get_entry("diffusive_jacobians",difje); CHKERRQ(difje==0);
    read_double_array(djacv,difje); 
    djacvp=djacv.begin();
    D_jac_l.set(0.);

    //o Scale the SUPG upwind term. 
    EGETOPTDEF_ND(elemset,double,tau_fac,1.);

    na = ajacv.size();
    nd = djacv.size();

    if (na==ndim*ndof && nd==ndof) {
      // An advection velocity and a diffusivity for each field 
      u.set(ajacv.begin());
      ret_options &= !SCALAR_TAU; // tell the advective element routine
				// that we are returning a non-scalar tau
      A_jac_l.set(0.);
      D_jac_l.set(0.);
      for (int k=1; k<=ndof; k++) {
	double alpha=djacv[k-1];
	for (int j=1; j<=ndim; j++) {
	  A_jac_l.setel(u.get(k,j),j,k,k);
	  D_jac_l.setel(alpha,j,j,k,k);
	}
      }

    } else if (na==ndim && (nd==1 || nd==ndof)) {
      // An advection velocity for all fields is entered
      assert(0);
#if 0 // Por ahora nos concentramos en el caso na==ndim*ndof
      vel = 0.;
      FastMat2 u(ndim),Uintri(ndim);
      Uintri.prod(iJaco,u,1,-1,-1);
      u.set(ajacv.begin());
      for (int j=1; j<=ndim; j++) {
	A_jac_l.ir(1,j).eye(ajacv[j-1]);
	vel += SQ(ajacv[j-1]);
      }
      vel=sqrt(vel);

      //        for (int k=1; k<=ndof; k++) {
      //  	double tau_a = tau_fac * h_supg/(2.* lam_max);
      //  	tau_supg.setel(
      //  		       }


      if (nd==1) {
	ret_options |= SCALAR_TAU; // tell the advective element routine
				// that we are returning a scalar tau
      } else {
	ret_options &= !SCALAR_TAU; // tell the advective element routine
				// that we are returning a scalar tau
      }
	
      A_jac_l.rs();

      h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
#endif

    } else if (na==ndim*ndof*ndof) {
      assert(0);
      // The full jacobian for each dimension
      A_jac_l.set(ajacv.begin());
      double vel=0.;
      for (int k=1; k<=ndim; k++) {
	A_jac_l.ir(1,k);
	AA.sum_abs(A_jac_l,1,-1);
	double a = AA.max_all();
	vel += a*a;
      }
      // I put a scalar tau by the moment
      vel = sqrt(vel);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Not a valid number of elements  while entering\n"
		  "the advective jacobians. Entered na=%d elements,\n"
		  "and ndim: %d, ndof %d\n"
		  "Allowed possibilities are: na=ndim, ndim*ndof, ndim*ndof*ndof\n",
		  na,ndim,ndof);
      PetscFinalize();
      exit(0);
    }
  }
  FastMat2::leave();

  A_jac.set(A_jac_l);
  flux.prod(A_jac,U,2,1,-1,-1);

  if (options & COMP_UPWIND) {

    D_jac.set(D_jac_l);
    fluxd.prod(D_jac,grad_U,2,-1,1,-2,-1,-2);

    // A_grad_U es ndof x 1
    A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);

    lam_max = 0.;
    tau_supg.set(0.);
    for (int k=1; k<=ndof; k++) {
      
      assert(na==ndim*ndof && nd==ndof);
      u.ir(1,k);
      double alpha=djacvp[k-1];
      Uintri.prod(iJaco,u,1,-1,-1);
      double vel = sqrt(u.sum_square_all());
      u.rs();
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
    delta_sc = 0.;
  }
  
  if (options & COMP_SOURCE) {
    G_source.set(0.);		// Only null source term allowed
				// right now!!
  }

}

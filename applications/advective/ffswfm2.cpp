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
#include "../../src/util2.h"

#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "flux_fun_shallowfm2" 
int flux_fun_shallowfm2(FLUX_FUN_ARGS_FM2) {

  static double ajacx[9],ajacy[9];
  int ierr;

  if (ndim!=2) {
    PFEMERRQ("Stop with \"mojito\" drinking!! 2D Only\n");
  }

  static int flag=0,ndof;
  static double g, gravity, tau_fac,shock_capturing_threshold;
  double tau_a, tau_delta, gU, A01v[9];
  static int shock_capturing;

  // Load properties only once.
  if (start_chunk) {
    start_chunk = 0;

    ndof = U.dim(1);
    //o Acceleration of gravity
    SGETOPTDEF_ND(double,gravity,1.);
    g=gravity;
    //o Scale the SUPG upwind term. 
    SGETOPTDEF_ND(double,tau_fac,1.);
    //o Add shock-capturing term.
    SGETOPTDEF_ND(int,shock_capturing,0);
    //o Add shock-capturing term if relative variation of variables
    // inside the element exceeds this.
    SGETOPTDEF_ND(double,shock_capturing_threshold,0.1);

    ajacx[1] = 0.;
    ajacx[6] = 1.;
    ajacx[7] = 0.;
    ajacx[8] = 0.;

    ajacy[3] = 0.;
    ajacy[6] = 0.;
    ajacy[7] = 1.;
    ajacy[8] = 0.;
  }

  static FMatrix u(ndim),flux_mass(ndim),UU,flux_mom,
    Uintri,vref(ndof),A01(ndof,ndof),tmp1,tmp2,tmp3,tmp4,grad_U_psi;
  UU.rs().set(U);

  double h = UU.get(ndim+1);
  double HH = H.get(1);
  flux_mass.set(UU.is(1,1,ndim));
  UU.rs();
  u.set(flux_mass).scale(1./h);

  double u2 = u.sum_square_all();

  double ux,uy;
  ux=u.get(1);
  uy=u.get(2);

  ajacx[0] = 2*ux;
  ajacx[2] = -ux*ux+g*h;
  ajacx[3] = uy;
  ajacx[4] = ux;
  ajacx[5] = -ux*uy;

  A_jac.ir(1,1).set(ajacx);

  ajacy[0] = uy;
  ajacy[1] = ux;
  ajacy[2] = -ux*uy;
  ajacy[4] = 2*uy;
  ajacy[5] = -uy*uy + g*h;

  A_jac.ir(1,2).set(ajacy).rs();

  flux_mom.prod(u,u,1,2).scale(h);

  // double h_term = 0.5*g*(h*h-HH*HH); !! ERROR
  double h_term = 0.5*g*h*h; 
  for (int jdim=1; jdim<=ndim; jdim++) {
    flux_mom.addel(h_term,jdim,jdim);
  }
  flux.rs().is(1,1,ndim).set(flux_mom);
  flux.rs().ir(1,ndof).set(flux_mass).rs();

  if (options & COMP_UPWIND) {

    // A_grad_U es ndof x 1
    A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);

    Uintri.prod(iJaco,u,1,-1,-1);
    double vel = sqrt(u2);
    double h_supg;

    FastMat2::branch();
    if (vel>1e-10) {
      FastMat2::choose(0);
      h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
    } else {
      // fixme:= This is for quads and hexas only?
      FastMat2::choose(1);
      double vol = 0.125/iJaco.det();
      h_supg = pow(vol,1./ndim);
    }
    FastMat2::leave();

    lam_max = fabs(sqrt(h*g)+vel);

    ret_options |= SCALAR_TAU;
    tau_a = tau_fac * h_supg/(2.* lam_max);

    double vmax = -1;

    vref.is(1,1,ndim).set(vel*h).rs();
    double pp = UU.get(ndof);
    vref.setel(pp,ndof);

    for (int jdof=1; jdof<=ndof; jdof++) {
           double vaux = vref.get(jdof);
	   vaux = vaux*vaux;
	   gU = grad_U.rs().ir(2,jdof).sum_square_all();
	   vaux = gU/vaux;
	   vmax = (vmax > vaux ? vmax : vaux);
    }
    grad_U.rs();
    vmax = sqrt(vmax);

    // Shock Capturing term. If not used return tau_supg as usual and
    // delta_sc=0.
    tau_delta = 0; delta_sc=0.;

    FastMat2::branch();
    if ( shock_capturing && (vmax > shock_capturing_threshold /h_supg) ) {
      FastMat2::choose(0);
      // calculo del tensor metrico de Riemann (A0) para transformar de variables 

      A01v[0]= 1.;
      A01v[1]= 0. ;
      A01v[2]= -ux;

      A01v[3]= 0.;
      A01v[4]= 1.;
      A01v[5]= -uy;

      A01v[6]= -ux;
      A01v[7]= -uy;
      A01v[8]= g*h+ux*ux+uy*uy;

      A01.set(A01v).scale(1./(g*h));

      // calculo del delta shock capturing delta_sc
      double vaux_num,vaux_den;
      tmp1.prod(A01,A_grad_U,1,-1,-1);
      tmp2.prod(A_grad_U,tmp1,-1,-1);
      vaux_num = double(tmp2);

      grad_U_psi.prod(iJaco,grad_U,1,-1,-1,2);
      tmp3.prod(A01,grad_U_psi,1,-1,2,-1);
      tmp4.prod(grad_U_psi,tmp3,-1,-2,-2,-1);
      vaux_den = double(tmp4);

      delta_sc = sqrt(vaux_num / vaux_den);
      tau_delta = delta_sc/(lam_max*lam_max);
      
    } // if (shock_capturing ...

    FastMat2::leave();

    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
    ret_options |= SCALAR_TAU;
    tau_supg.setel(tau_supg_d,1,1);

  } 

  if (options & COMP_SOURCE) {
    G_source.set(0.);
    G_source.is(1,1,ndim).set(grad_H.ir(2,1)).scale(-g*h).rs();
    grad_H.rs();
  }

}

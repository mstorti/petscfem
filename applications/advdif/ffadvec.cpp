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
#include "newmat.h"
#include <vector>

#include "../../src/fem.h"
#include "../../src/texthash.h"
#include "../../src/getprop.h"

#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int flux_fun_advec()"
int flux_fun_advec(FLUX_FUN_ARGS) {

  int ierr;

  if (ndim!=2) {
    PFEMERRQ("Stop with \"mojito\" drinking!! 2D Only\n");
  }

  int ndof = U.Ncols();
  double vell;

  static int flag=0;
  static double vel=0, tau_fac=1.;
  static Matrix u;
  static ColumnVector uu;
  int kmax;

  if (flag==0) {
    uu.ReSize(ndim);
    flag=1;
    u.ReSize(ndim,ndof);
    u << 1. << 0 << 0.5 << 0. << -1. << 0.;
    for (int k=1; k<=ndof; k++) {
      vell = (u.Column(k)).SumSquare();
      if (k==1 || vell> vel) {
	vel = vell;
	kmax = k;
      }
    }
    for (int jd=1; jd<=ndim; jd++) {
      uu(jd) = u(jd,kmax);
    }
  }    
  
  for (int jd=1; jd<=ndim; jd++) {
    AJAC(jd)= 0;
    for (int k=1; k<=ndof; k++) {
      AJAC(jd)(k,k) = u(jd,k);
    }
  }

  // A_grad_U es ndof x 1
  A_grad_U = 0;
  for (int jd=1; jd<=ndim; jd++) {
    A_grad_U += AJAC(jd) * grad_U.Row(jd).t();
  }

  if ( options & COMP_UPWIND ) {
    Matrix Uintri = iJaco * uu;
    double h_supg = 2.*vel/sqrt(Uintri.SumSquare());
    
    tau_supg = tau_fac * h_supg/(2.* vel);
    lam_max = vel;
  }
	
  if (options & COMP_SOURCE) G_source = 0.;

}



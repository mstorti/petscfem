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
#define __FUNC__ "flux_fun_euler" 
int flux_fun_euler(FLUX_FUN_ARGS) {

  // Convert to FastMat
#define CONV_TO_FM(M) FastMat M##_FM; NM2FM(M##_FM,M)
  CONV_TO_FM(U);
  CONV_TO_FM(iJaco);
  CONV_TO_FM(H);
  CONV_TO_FM(grad_H);
  CONV_TO_FM(grad_U);
  CONV_TO_FM(nor);
  CONV_TO_FM(lambda);
  CONV_TO_FM(Vr);
  CONV_TO_FM(Vr_inv);
#undef CONV_TO_FM

  vector<FastMat *> A_jac_FM;
  for (int jd=1; jd<=ndim; jd++)
    A_jac_FM.push_back(new FastMat(ndof,ndof));

  flux_fun_euler_FM(FLUX_FUN_CALL_ARGS_FM);

#define CONV_TO_NM(M) M << M##_FM.store  
  // Convert back to Newmat
  CONV_TO_NM(flux);
  CONV_TO_NM(A_grad_U);
  CONV_TO_NM(G_source);
  CONV_TO_NM(tau_supg);
  CONV_TO_NM(lambda);
  CONV_TO_NM(Vr);
  CONV_TO_NM(Vr_inv);
#undef CONV_TO_NM

  vector<FastMat *> A_jac_FM;
  for (int jd=1; jd<=ndim; jd++) {
    *(A_jac[jd-1]) << A_jac_FM.store;
    delete A_jac_FM[jd-1];
  }
  
}

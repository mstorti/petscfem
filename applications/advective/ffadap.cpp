//__INSERT_LICENSE__
//$Id: ffadap.cpp,v 1.2 2001/04/01 01:34:41 mstorti Exp $

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

//__INSERT_LICENSE__
//$Id: ffadapfm2.cpp,v 1.5 2001/12/20 21:58:55 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "advective.h"

int flux_fun_euler_FM2(FLUX_FUN_ARGS_FM2);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "flux_fun_euler" 
int flux_fun_eulerfm2(FLUX_FUN_ARGS) {

#if 1
  static FastMatCacheList cache_list;
  if (start_chunk) {
    FastMat2::void_cache();
    FastMat2::activate_cache(&cache_list);
    start_chunk=0;
  }
  FastMat2::reset_cache();
#endif

  // Convert to FastMat
#define CONV_TO_FM2(M) static FastMat2 M##_FM2; M##_FM2.set(M)
  CONV_TO_FM2(U);
  CONV_TO_FM2(flux);
  CONV_TO_FM2(tau_supg);
  CONV_TO_FM2(iJaco);
  CONV_TO_FM2(H);
  CONV_TO_FM2(grad_H);
  CONV_TO_FM2(grad_U);
  CONV_TO_FM2(nor);
  CONV_TO_FM2(lambda);
  CONV_TO_FM2(Vr);
  CONV_TO_FM2(Vr_inv);
  CONV_TO_FM2(G_source);
#undef CONV_TO_FM2

  int ndof = U.Ncols();
  static FastMat2 A_grad_U_FM2;
  static FastMat2 A_jac_FM2(3,ndim,ndof,ndof);

  flux_fun_euler_FM2(FLUX_FUN_CALL_ARGS_FM2);

  //#define CONV_TO_NM(M) FM2_2_NM(M,M##_FM2)
#define CONV_TO_NM(M) M##_FM2.export(M)
  // Convert back to Newmat
  CONV_TO_NM(flux);
  CONV_TO_NM(A_grad_U);
  CONV_TO_NM(G_source);
  CONV_TO_NM(tau_supg);
  CONV_TO_NM(lambda);
  CONV_TO_NM(Vr);
  CONV_TO_NM(Vr_inv);    
#undef CONV_TO_NM

  // static FastMat2 A_jac_j(2,ndof,ndof);
  for (int jd=1; jd<=ndim; jd++) {
    A_jac_FM2.rs().ir(1,jd).export(*(A_jac[jd-1]));
//      A_jac_j.set(A_jac_FM2.rs().ir(1,jd));
//      FM2_2_NM(*(A_jac[jd-1]),A_jac_j);
  }
  
}

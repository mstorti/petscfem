/*__INSERT_LICENSE__*/
//$Id: ffadvec.cpp,v 1.3 2001/05/30 03:58:38 mstorti Exp $

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



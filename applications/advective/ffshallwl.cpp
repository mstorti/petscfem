//__INSERT_LICENSE__
//$Id: ffshallwl.cpp,v 1.4.10.1 2001/12/21 00:13:32 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include "newmat.h"
#include <vector>
#include <cassert>

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include "advective.h"
#include "util2.h"

double inv_abs(double r) {
  return 1./fabs(r);
}

int flux_fun(const RowVector &U,int ndim,const Matrix &iJaco,
	     Matrix &H, Matrix &grad_H,
	     Matrix &flux,vector<Matrix *> A_jac,Matrix &G_source,Matrix &tau_supg,
	     TextHashTable *thash,double *propel,void *user_data,
	     int options) {

#define NON_ISOTROPIC_TAU

  int ierr;

  if (ndim!=2) {
    PFEMERRQ("Stop with \"mojito\" drinking!! 2D Only\n");
  }

  int ndof = U.Ncols();

  static int flag=0;
  static double g, tau_fac=1.;

  // Load properties only once.
  if (flag==0) {
    flag = 1;
    ierr = get_double(thash,"gravity",&g); CHKERRA(ierr);
  }

  double h,ux,uy;
  double HH = H.AsScalar();

  h = 1.;
  ux=0.3;
  uy=0;

  AJAC(1) << 2*ux << 0 << -ux*ux+g*h 
     << uy << ux << -ux*uy 
     << 1  << 0 << 0 ;

#if 0
  AJAC(2) <<   uy << ux << -ux*uy
     << 0  << 2*uy << -uy*uy + g*h
     << 0 << 1  << 0 ;
#else
  // fixme:= para hacerlo unidimensional
  AJAC(2) = 0;
#endif

  if (options & COMP_UPWIND) {

    ColumnVector u(2);
    u(1)=ux;
    u(2)=uy;
    Matrix Uintri = iJaco * u;
    double vel = sqrt(u.SumSquare());
    double h_supg = 2.*vel/sqrt(Uintri.SumSquare());
    double lam_max = vel + sqrt(g*h);

#ifdef NON_ISOTROPIC_TAU
    tau_supg.ReSize(ndof,ndof);
    mat_function(AJAC(1),tau_supg,&inv_abs);
    tau_supg = (tau_fac * h_supg / 2.) * tau_supg;
#else
    tau_supg.ReSize(1,1);
    tau_supg = tau_fac * h_supg/(2.* lam_max);
#endif
    
	
  } 

  if (options & COMP_SOURCE) {

    G_source = 0.;
    // fixme:= Esto esta en discusion. 
    G_source.Rows(1,ndim) =  -(g * (h+HH)) * grad_H;
    // G_source.Rows(1,ndim) =  (g * (h-HH)) * grad_H;
  }

}

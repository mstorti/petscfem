//__INSERT_LICENSE__
//$Id: ffshallw.cpp,v 1.2 2001/04/01 01:34:41 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include "../../src/fem.h"
#include "../../src/texthash.h"
#include "../../src/getprop.h"
#include "../../src/utils.h"

#include "advective.h"

double inv_abs(double r,void * user_data) {
  double lam_cut = *(double *) user_data;
  double la;
  la = (fabs(r) < lam_cut ? lam_cut : fabs(r));
  return 1./la;
}

ColumnVector tau_fun(ColumnVector v,void * user_data) {
  double addvisc = *(double *) user_data;
  double md = 0;
  int n = v.Nrows();
  ColumnVector fv(n);

  for (int k=1; k<=v.Nrows(); k++) {
    fv(k) = fabs(v(k));
    md = (fv(k)>md ? fv(k) : md );
  }

  for (int k=1; k<=fv.Nrows(); k++) {
    fv(k) = (fv(k)+ addvisc *(md-fv(k)))/(fv(k)*fv(k));
  }
  return fv;
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "flux_fun" 
int flux_fun_shallow(FLUX_FUN_ARGS) {

  ret_options |= SCALAR_TAU;
  delta_sc=0.; // no shock_capturing for shallow water!!
  int ierr;

  if (ndim!=2) {
    PFEMERRQ("Stop with \"mojito\" drinking!! 2D Only\n");
  }

  int ndof = U.Ncols();

  static int flag=0;
  static double g, tau_fac=1., addvisc;

  // Load properties only once.
  if (flag==0) {
    flag = 1;
    ierr = get_double(thash,"gravity",&g); CHKERRA(ierr);
    // ierr = get_double(thash,"addvisc",&addvisc); CHKERRA(ierr);
  }

  ColumnVector u;
  RowVector flux_mass;
  double h = U(ndim+1);
  double HH = H.AsScalar();
  flux_mass = U.Columns(1,ndim);
  u = flux_mass.t()/h;
  double ux,uy;
  ux=u(1);
  uy=u(2);

  AJAC(1) << 2*ux << 0 << -ux*ux+g*h 
     << uy << ux << -ux*uy 
     << 1  << 0 << 0 ;

  //#define UPWIND_1D

#ifdef UPWIND_1D
  AJAC(2) = 0;
#else
  AJAC(2) <<   uy << ux << -ux*uy
     << 0  << 2*uy << -uy*uy + g*h
     << 0 << 1  << 0 ;
#endif

  Matrix flux_mom = h * u * u.t();

  // double h_term = 0.5*g*(h*h-HH*HH); !! ERROR
  double h_term = 0.5*g*h*h; 
  for (int jdim=1; jdim<=ndim; jdim++) {
    flux_mom(jdim,jdim) += h_term;
  }
  flux.SubMatrix(1,ndim,1,ndim) = flux_mom;
  flux.Row(ndof) = flux_mass;

  if (options & COMP_UPWIND) {

    // A_grad_U es ndof x 1
    A_grad_U = 0;
    for (int jd=1; jd<=ndim; jd++) {
      A_grad_U += AJAC(jd) * grad_U.Row(jd).t();
    }

    Matrix Uintri = iJaco * u;

    double vel = sqrt(u.SumSquare());
    double h_supg;

    if (vel>1e-10) {
      h_supg = 2.*vel/sqrt(Uintri.SumSquare());
    } else {
      double vol = 0.125/mydet(iJaco);
      h_supg = pow(vol,1./ndim);
    }


#ifdef UPWIND_1D
    tau_supg.ReSize(ndof,ndof);
    lam_max = vel+ sqrt(h*g);
    // lam_cut = .3*lam_max;
    //    lam_cut = 0.;
    mat_function(AJAC(1),tau_supg,&tau_fun,(void *)&addvisc);
    tau_supg = (tau_fac * h_supg / 2.) * tau_supg;
    
#else
    //guarda esto no se como esta
    lam_max = fabs(sqrt(h*g)+vel);
    tau_supg(1,1) = tau_fac * h_supg/(2.* lam_max);
#endif
    
	
  } 

  if (options & COMP_SOURCE) {

    G_source = 0.;
    // fixme:= Esto esta en discusion. 
    G_source.Rows(1,ndim) =  -g * h * grad_H;
    // Si H esta definido como positivo hacia abajo. 
    // G_source.Rows(1,ndim) =  (g * (h-HH)) * grad_H;
  }

}

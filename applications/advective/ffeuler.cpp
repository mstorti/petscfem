//__INSERT_LICENSE__
//$Id: ffeuler.cpp,v 1.2 2001/04/01 01:34:41 mstorti Exp $

// THIS FILE APPARENTLY DOESN'T COMPILE WITH OPTIMIZATION.
// SO THAT I MODIFIED THE RULE IN THE MAKEFILE SO THAT IT
// IS ALWAYS COMPILED WITH `-O0'

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include "../../src/fem.h"
#include "../../src/texthash.h"
#include "../../src/getprop.h"
#include "../../src/util2.h"

#include "advective.h"

//#define USE_FASTMAT

#ifndef USE_FASTMAT
#define MATRIX_LIB Matrix
#else
#define MATRIX_LIB FastMat
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "flux_fun_euler" 
#ifndef USE_FASTMAT
int flux_fun_euler(FLUX_FUN_ARGS) {
#else
int flux_fun_euler_FM(FLUX_FUN_ARGS_FM) {
#endif

  int ierr;
  // fixme:= 
  // assert(0); // no se porque no compila, lo comento por ahora
  ret_options |= SCALAR_TAU;
#ifndef USE_FASTMAT
  int ndof = U.Ncols();
#else
  int ndof = U.n;
#endif

  assert(ndof==ndim+2);

  static int flag=0, shock_capturing;
  static double tau_fac, gamma;
  double tau_a, tau_delta;

  // Load properties only once.
  if (flag==0) {
    flag = 1;
    gamma=1.4;
    ierr = get_double(thash,"gamma",&gamma,1); CHKERRA(ierr);
    tau_fac=1.;
    ierr = get_double(thash,"tau_fac",&tau_fac,1); CHKERRA(ierr);
    shock_capturing=0;
    ierr = get_int(thash,"shock_capturing",&shock_capturing,1); CHKERRA(ierr);
  }
  MATRIX_LIB A01(ndof,ndof),grad_U_psi(ndim,ndof),
    Omega(3,3),VV3d(5,5),VV3di(5,5),VV(ndof,ndof),VVi(ndof,ndof),
    vaux,waux,Uintri;

#ifndef USE_FASTMAT
  ColumnVector u,v(ndim),u3d(3),kdir(ndim),kdir3d(3);
  RowVector V(ndof),vref(ndof);
  DiagonalMatrix I(ndim);
  I=1;
#else
  FastMat u,u_t(1,ndim),v(ndim,1),u3d(3,1),V(1,ndof),kdir(ndim,1),kdir3d(3,1),I,
    vref(1,ndof),kdir3d_t(1,3),u3d_t(1,3);
  I.eye(ndim);
  FastMat tmp1,tmp2,tmp3,tmp5,tmp6,tmp7,tmp13,
    tmp15,tmp17,tmp8,tmp10,tmp18,tmp19,tmp20,tmp21;
#endif

  // fixme:= vaux esta usado de varias formas segun el scope
  double p,energy;

#ifndef USE_FASTMAT
  double rho = U(1);
  energy = U(ndim+2)/rho; // total energy
  u = (U.Columns(2,ndim+1)/rho).t();
#else
  double rho = U.get(1,1);
  energy = U.get(1,ndim+2)/rho; // total energy
  u_t.columns(U,2,ndim+1);
  u_t.scale(1/rho);
  u.transpose(u_t);
#endif

  double g1=gamma-1;
#ifndef USE_FASTMAT
  double u2 = u.SumSquare(); 
#else
  double u2;
  u.sum_square(u2); 
#endif

  double i = energy - 0.5 * u2; // internal energy
  p = g1*rho*i; // pressure
  double cc = sqrt(gamma*p/rho);

#ifndef USE_FASTMAT
  Matrix Z(1,1), ZZ(1,1), ZZZ(1,1);
  Z=0.;
#endif

  for (int jdim=1; jdim<=ndim; jdim++) {

#ifndef USE_FASTMAT
    double uk = u(jdim);
    ZZ = uk*(g1*u2-gamma*energy);
    ZZZ = gamma*uk;
    kdir=0;
    kdir(jdim)= 1;
    AJAC(jdim) = (Z | kdir.t() | Z) &
      (-uk*u + 0.5*g1*u2*kdir | u*kdir.t()+uk*I-g1*kdir*u.t() | g1*kdir) &
      (ZZ | (p/rho+energy)*kdir.t()-g1*uk*u.t() | ZZZ);
    flux.Column(jdim) = AJAC(jdim) * U.t();
#else
    double uk = u.get(jdim,1);
    kdir.set(0.);
    kdir.set(jdim,1,1.);

    // the corners
    AJAC(jdim).set(1,1,0.);
    AJAC(jdim).set(1,ndim+2,0.);
    AJAC(jdim).set(ndof,1,uk*(g1*u2-gamma*energy));
    AJAC(jdim).set(ndof,ndof,gamma*uk);
    
    // vector like submatrices
    for (int jd=1; jd<=ndim; jd++) {
      AJAC(jdim).set(1,jd+1,kdir.get(jd,1));
      AJAC(jdim).set(jd+1,1,-uk*u.get(jd,1) + 0.5*g1*u2*kdir.get(jd,1));
      AJAC(jdim).set(jd+1,ndof,g1*kdir.get(jd,1));
      AJAC(jdim).set(ndof,jd+1,(p/rho+energy)*kdir.get(jd,1)-g1*uk*u.get(jd,1));

      // the core
      for (int kd=1; kd<=ndim; kd++) {
	AJAC(jdim).set(jd+1,kd+1,
		       u.get(jd,1)*kdir.get(kd,1)+
		       uk*I.get(jd,kd)-g1*kdir.get(jd,1)*u.get(kd,1));
      }
    }
#endif    
  }

  if ( options & COMP_UPWIND ) {

    // A_grad_U es ndof x 1

#ifndef USE_FASTMAT
    A_grad_U = 0;
#else
    A_grad_U.set(0.);
#endif

    for (int jd=1; jd<=ndim; jd++) {
#ifndef USE_FASTMAT
      A_grad_U += AJAC(jd) * grad_U.Row(jd).t();
#else
      tmp1.row(grad_U,jd);
      tmp2.transpose(tmp1);
      FMp(tmp3,AJAC(jd),tmp2);
      FMaxpy(A_grad_U,1.,tmp3);
#endif
    }

#ifndef USE_FASTMAT
    Uintri = iJaco * u;
#else
    FMp(Uintri,iJaco,u);
#endif
    double vel = sqrt(u2);

#ifndef USE_FASTMAT
    double h_supg = 2.*vel/sqrt(Uintri.SumSquare());
#else
    double uintri2;
    Uintri.sum_square(uintri2);
    double h_supg = 2.*vel/sqrt(uintri2);
#endif

    double c = sqrt(gamma*p/rho);
#if 0
#ifndef USE_FASTMAT
    tau_supg.ReSize(1,1);
#else
    tau_supg.reshape(1,1);
#endif
#endif

    lam_max = c+vel;
    // double lam_max = c;
    tau_a = tau_fac * h_supg/(2.* lam_max);

    double vmax = -1;

#ifndef USE_FASTMAT
    vref(1) = U(1);
    vref.Columns(2,ndim+1) = rho*vel;
    vref(ndim+2)=U(ndof);
#else
    vref.set(1,1,U.get(1,1));
    vref.set(2,ndim+1,rho*vel);
    vref.set(ndof,ndof,U.get(1,ndof));
#endif

    for (int jdof=1; jdof<=ndof; jdof++) {

#ifndef USE_FASTMAT
           double vaux = pow(vref(jdof),2);
	   vaux = (grad_U.Column(jdof)).SumSquare()/vaux;
#else
	   double tmp4,aux6;
	   tmp4 = vref.get(1,jdof);
           double vaux = tmp4*tmp4;
	   tmp5.columns(grad_U,jdof,jdof);
	   tmp5.sum_square(aux6);
	   vaux = aux6/vaux;
#endif
	   vmax = (vmax > vaux ? vmax : vaux);
    }
    vmax = sqrt(vmax);

    // Shock Capturing term. If not used return tau_supg as usual and
    // delta_sc=0.
    tau_delta = 0; delta_sc=0.;
    if ( shock_capturing && (vmax > 0.1/h_supg) ) {
      // calculo del tensor metrico de Riemann (A0) para transformar de variables 
      // conservativas a variables de entropia (ver Mallet o Shakib)
#ifndef USE_FASTMAT
      V(1) = 0; 
#else
      V.set(1,1,0.);
#endif
      for (int jdim=1; jdim<=ndim; jdim++) {
#ifndef USE_FASTMAT
	V(jdim+1) = U(jdim+1)/(rho*i);
#else
	V.set(1,jdim+1,U.get(1,jdim+1)/(rho*i));
#endif
      }

#ifndef USE_FASTMAT
      V(ndof) = -U(1)/(rho*i);
#else
      V.set(1,ndof,-U.get(1,1)/(rho*i));
#endif

#ifndef USE_FASTMAT
      v = V.Columns(2,ndim+1).t(); 
      double k1  = 0.5*SumSquare(v)/V(ndof);
      ZZ = k1*k1+gamma;
      ZZZ = (k1+1)*V(ndof);
      Z = V(ndof)*V(ndof);
      A01 = (ZZ | k1*v.t() | ZZZ ) &
	(k1*v | v*v.t()-V(ndof)*I | V(ndof)*v ) &
	(ZZZ | V(ndof)*v.t() | Z );

      A01 = -A01/(rho*i*V(ndof));
#else
      double Vndof = V.get(1,ndof);
      v.columns(V,2,ndim+1);
      double aux7,k1;
      v.sum_square(aux7);
      k1=0.5*aux7/Vndof;
      A01.set(1,1,k1*k1+gamma);
      double aux8=(k1+1)*Vndof;
      A01.set(1,ndof,aux8);
      A01.set(ndof,1,aux8);
      A01.set(ndof,ndof,Vndof*Vndof);
      for (int jd=1; jd<=ndim;jd++) {

	A01.set(1,jd+1,k1*v.get(jd,1));
	A01.set(jd+1,1,k1*v.get(jd,1));

	A01.set(ndof,jd+1,Vndof*v.get(jd,1));
	A01.set(jd+1,ndof,Vndof*v.get(jd,1));

	for (int kd=1; kd<=ndim; kd++) 
	  A01.set(jd,kd,v.get(jd,1)*v.get(kd,1));
	A01.add(jd,jd,-Vndof);
      }
#endif

    // calculo del delta shock capturing delta_sc
      double vaux_num,vaux_den;
#ifndef USE_FASTMAT
      vaux_num = (A_grad_U.t()*A01*A_grad_U).AsScalar();
      grad_U_psi = iJaco * grad_U;
      vaux_den = (grad_U_psi*A01*grad_U_psi.t()).Trace();
#else
      FMp(tmp2,A01,A_grad_U);
      tmp1.transpose(A_grad_U);
      tmp1.trace_of_product(tmp2,vaux_num);

// tmp6 es de ndim x ndof y 
// tmp7 y tmp8 es de ndof x ndim

      FMp(tmp6,iJaco,grad_U);
      tmp7.transpose(tmp6);
      FMp(tmp8,A01,tmp7);
      tmp6.trace_of_product(tmp8,vaux_den);
#endif
      delta_sc = sqrt(vaux_num / vaux_den);
      tau_delta = delta_sc/(lam_max*lam_max);

    }
    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
#ifndef USE_FASTMAT
    tau_supg(1,1) = tau_supg_d;
#else
    tau_supg.set(1,1,tau_supg_d);
#endif
  } 

  if (options & COMP_SOURCE) {

#ifndef USE_FASTMAT
    G_source = 0.;
#else
    G_source.set(0.);
#endif
  }

  if (options & COMP_EIGENV) {

#ifndef USE_FASTMAT
    kdir3d=0;
    kdir3d.Rows(1,ndim)=nor;
#else
    kdir3d.set(0.);
    kdir3d.set(1,ndim,1,1,nor);
    kdir3d_t.transpose(kdir3d);
#endif

#ifndef USE_FASTMAT
    u3d = 0;
    Z = 0.;
    u3d.Rows(1,ndim) = u;
#else
    u3d.set(0.);
    u.rows(u3d,1,ndim);
    u3d_t.transpose(u3d);
#endif

#ifndef USE_FASTMAT
    Omega << 0. << -kdir3d(3) << kdir3d(2) <<
      kdir3d(3) << 0. << -kdir3d(1) <<
      -kdir3d(2) << kdir3d(1) << 0.;
#else
    Omega.set(0.);
    Omega.set(1,2,-kdir3d.get(3,1));
    Omega.set(1,3,kdir3d.get(2,1));
    Omega.set(2,1,kdir3d.get(3,1));
    Omega.set(2,3,-kdir3d.get(1,1));
    Omega.set(3,1,-kdir3d.get(2,1));
    Omega.set(3,2,kdir3d.get(1,1));
#endif

    double uk,z;
#ifndef USE_FASTMAT
    uk = (u3d.t()*kdir3d).AsScalar();
#else
    u3d_t.trace_of_product(kdir3d,uk);
#endif
    z=0.5*sqrt(2.0)*rho/cc;

#ifndef USE_FASTMAT
    Z = z;
    VV3d = (kdir3d.t() | Z | Z) &
      (u3d*kdir3d.t()+rho*Omega | z*(u3d+cc*kdir3d) | z*(u3d-cc*kdir3d) ) &
      ( 0.5*u2*kdir3d.t()-rho*(Omega*u3d).t() | 
	Z*(0.5*u2+cc*cc/g1+cc*uk) | Z*(0.5*u2+cc*cc/g1-cc*uk) );
#else
    tmp2.reshape(3,1);
    FMp(tmp2,Omega,u3d);
    for (int jd=1; jd<=3; jd++) {
      VV3d.set(1,jd,kdir3d_t.get(1,jd));
      VV3d.set(5,jd,0.5*u2*kdir3d_t.get(1,jd)-rho*tmp2.get(jd,1));
      VV3d.set(jd+1,4,z*(u3d.get(jd,1)+cc*kdir3d.get(jd,1)));
      VV3d.set(jd+1,5,z*(u3d.get(jd,1)-cc*kdir3d.get(jd,1)));
      for (int id=1; id<=3; id++) {
        VV3d.set(jd+1,id,
              u3d.get(jd,1)*kdir3d.get(id,1)+rho*Omega.get(jd,id));
      }
    }
    VV3d.set(1,4,z);
    VV3d.set(1,5,z);
    VV3d.set(5,4,z*(0.5*u2+cc*cc/g1+cc*uk));
    VV3d.set(5,5,z*(0.5*u2+cc*cc/g1-cc*uk));
#endif

    z = 0.5*sqrt(2.0)/rho/cc;
#ifndef USE_FASTMAT
    Z = z;
    VV3di = ( (1-0.5*g1*u2/cc/cc)*kdir3d+(1./rho)*(Omega*u3d) | 
	      g1/cc/cc*kdir3d*u3d.t()-(1./rho)*Omega | -g1/cc/cc*kdir3d ) &
      (Z*(-cc*uk+0.5*g1*u2) | z*( cc*kdir3d.t()-g1*u3d.t()) | Z*g1 ) &
      (Z*( cc*uk+0.5*g1*u2) | z*(-cc*kdir3d.t()-g1*u3d.t()) | Z*g1 );
#else
    for (int jd=1; jd<=3; jd++) {
      VV3di.set(jd,1,
          (1-0.5*g1*u2/cc/cc)*kdir3d.get(jd,1)+(1./rho)*tmp2.get(jd,1));
      VV3di.set(jd,5,-g1/cc/cc*kdir3d.get(jd,1));
      VV3di.set(4,jd+1,z*( cc*kdir3d.get(jd,1)-g1*u3d.get(jd,1)));
      VV3di.set(5,jd+1,z*(-cc*kdir3d.get(jd,1)-g1*u3d.get(jd,1)));
      for (int id=1; id<=3; id++) {
        VV3di.set(jd,id+1,
          g1/cc/cc*kdir3d.get(jd,1)*u3d.get(id,1)-(1./rho)*Omega.get(jd,id));
      }
    }
    VV3di.set(4,1,z*(-cc*uk+0.5*g1*u2));
    VV3di.set(4,5,z*g1);
    VV3di.set(5,1,z*(cc*uk+0.5*g1*u2));
    VV3di.set(5,5,z*g1);
#endif


#ifndef USE_FASTMAT
    lambda(1,1) = uk; lambda(2,1) = uk;
#else
    lambda.set(1,1,uk); lambda.set(2,1,uk);
#endif

    if (ndim==2) {
#ifndef USE_FASTMAT
      lambda(3,1) = uk+cc; lambda(4,1) = uk-cc;
#else
      lambda.set(3,1,uk+cc); lambda.set(4,1,uk-cc);
#endif

#ifndef USE_FASTMAT
      vaux = (VV3d.Column(1)+VV3d.Column(3));
      Vr.Column(1) = (vaux.Rows(1,3) & vaux.Row(5)) ;
      vaux = (VV3d.Column(2)-VV3d.Column(3));	
      Vr.Column(2) = (vaux.Rows(1,3) & vaux.Row(5)) ;
#else
      tmp3.reshape(5,1);
      tmp5.reshape(5,1);
      tmp6.reshape(5,1);

      tmp3.column(VV3d,1);
      tmp5.column(VV3d,3);
      FMaxpy(tmp3,1.,tmp5);
      tmp10.rows(tmp3,1,3);
      Vr.set(1,3,1,1,tmp10);
      Vr.set(4,1,tmp3.get(5,1));

      tmp3.column(VV3d,2);
      FMaxpy(tmp3,-1.,tmp5);
      tmp10.rows(tmp3,1,3);
      Vr.set(1,3,2,2,tmp10);
      Vr.set(4,2,tmp3.get(5,1));
#endif

#ifndef USE_FASTMAT
      vaux.ReSize(5,2);
      vaux = VV3d.Columns(4,5);
      Vr.Columns(3,4) = vaux.Rows(1,3) & vaux.Row(5);
#else
      tmp17.get(1,3,4,5,VV3d);
      Vr.set(1,3,3,4,tmp17);
      tmp18.get(5,5,4,5,VV3d);
      Vr.set(4,4,3,4,tmp18);
#endif

#ifndef USE_FASTMAT
      if (kdir3d(1)*kdir3d(2)!=0.) {
	waux = (VV3di.Row(1)/kdir3d(1)+VV3di.Row(3)*kdir3d(2)) / 
	  (kdir3d(1)+kdir3d(2));
	Vr_inv.Row(1) = (waux.Columns(1,3) | waux.Column(5)) ; 
	waux = (VV3di.Row(2)/kdir3d(2)-VV3di.Row(3)*kdir3d(1)) / 
	  (kdir3d(1)+kdir3d(2));
	Vr_inv.Row(2) = (waux.Columns(1,3) | waux.Column(5)) ; 
      }
      else if (kdir3d(1)!=0) {
	waux = VV3di.Row(1);
	Vr_inv.Row(1) = (waux.Columns(1,3) | waux.Column(5)) ; 
	waux = VV3di.Row(1)-VV3di.Row(3);
	Vr_inv.Row(2) = (waux.Columns(1,3) | waux.Column(5)) ; 
      }
      else if (kdir3d(2)!=0) {
	waux = VV3di.Row(2)+VV3di.Row(3);
	Vr_inv.Row(1) = (waux.Columns(1,3) | waux.Column(5)) ; 
	waux = VV3di.Row(2);
	Vr_inv.Row(2) = (waux.Columns(1,3) | waux.Column(5)) ; 
      }
#else
      if (kdir3d.get(1,1)*kdir3d.get(2,1)!=0.) {
        tmp13.row(VV3di,1);
        tmp15.row(VV3di,3);
        tmp13.scale(1/kdir3d.get(1,1));
        tmp15.scale(kdir3d.get(2,1));
        FMaxpy(tmp13,1.,tmp15);
        tmp13.scale(1/(kdir3d.get(1,1)+kdir3d.get(2,1)));
	tmp19.get(1,1,1,3,tmp13);
        Vr_inv.set(1,1,1,3,tmp19);
        Vr_inv.set(1,1,4,4,tmp13.get(1,5));

        tmp13.row(VV3di,2);
        tmp15.row(VV3di,3);
        tmp13.scale(1/kdir3d.get(2,1));
        tmp15.scale(kdir3d.get(1,1));
        FMaxpy(tmp13,-1.,tmp15);
        tmp13.scale(1/(kdir3d.get(1,1)+kdir3d.get(2,1)));
	tmp19.get(1,1,1,3,tmp13);
        Vr_inv.set(2,2,1,3,tmp19);
        Vr_inv.set(2,2,4,4,tmp13.get(1,5));
      }
      else if (kdir3d.get(1,1)!=0) {
        tmp13.row(VV3di,1);
	tmp19.get(1,1,1,3,VV3di);
        Vr_inv.set(1,1,1,3,tmp19);
        Vr_inv.set(1,1,4,4,VV3di.get(1,5));
        tmp15.row(VV3di,3);
        FMaxpy(tmp13,-1.,tmp15);
	tmp19.get(1,1,1,3,tmp13);
        Vr_inv.set(2,2,1,3,tmp19);
        Vr_inv.set(2,2,4,4,tmp13.get(1,5));
      }
      else if (kdir3d.get(2,1)!=0) {
        tmp13.row(VV3di,2);
	tmp19.get(1,1,1,3,tmp13);
        Vr_inv.set(2,2,1,3,tmp19);
        Vr_inv.set(2,2,4,4,tmp13.get(1,5));
        tmp15.row(VV3di,3);
        FMaxpy(tmp13,1.,tmp15);
	tmp19.get(1,1,1,3,tmp13);
        Vr_inv.set(1,1,1,3,tmp19);
        Vr_inv.set(1,1,4,4,tmp13.get(1,5));
      }

#endif
#ifndef USE_FASTMAT
//      waux.ReSize(5,2);
      waux = VV3di.Rows(4,5);
      Vr_inv.Rows(3,4) = waux.Columns(1,3) | waux.Column(5);
#else
      tmp20.get(4,5,1,3,VV3di);
      Vr_inv.set(3,4,1,3,tmp20);
      tmp21.get(4,5,5,5,VV3di);
      Vr_inv.set(3,4,4,4,tmp21);
#endif

    }  else {
#ifndef USE_FASTMAT
      lambda(3,1) = uk; lambda(4,1) = uk+cc; lambda(5,1) = uk-cc;
      Vr = VV3d; Vr_inv = VV3di;
#else
      lambda.set(3,1,uk); lambda.set(4,1,uk+cc); lambda.set(5,1,uk-cc);
      Vr.set(VV3d); Vr_inv.set(VV3di);
#endif
    }
  }
}

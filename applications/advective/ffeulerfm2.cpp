//__INSERT_LICENSE__
//$Id: ffeulerfm2.cpp,v 1.7 2002/01/14 03:45:05 mstorti Exp $

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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "int flux_fun_euler()" 
int flux_fun_eulerfm2(FLUX_FUN_ARGS_FM2) {

  int ierr;
  int ndof = U.dim(1);

  assert(ndof==ndim+2);

  static int flag=0, shock_capturing;
  static double tau_fac, gamma,shock_capturing_threshold;
  double tau_a, tau_delta;

  // Load properties only once.
  if (start_chunk) {
    flag = 1;
    //o The specific heat ratio. 
    SGETOPTDEF_ND(double,gamma,1.4);
    //o Scale the SUPG upwind term. 
    SGETOPTDEF_ND(double,tau_fac,1.);
    //o Add shock-capturing term.
    SGETOPTDEF_ND(int,shock_capturing,0);
    //o Add shock-capturing term if relative variation of variables
    // inside the element exceeds this.
    SGETOPTDEF_ND(double,shock_capturing_threshold,0.1);
    CHKERRA(ierr);
    start_chunk=0;
  }
  static FastMat2 A01(2,ndof,ndof),grad_U_psi(2,ndim,ndof),
    Omega(2,3,3),VV3d(2,5,5),VV3di(2,5,5),VV(2,ndof,ndof),VVi(2,ndof,ndof),
    vaux,waux,Uintri,UU(1,ndof);

  static FastMat2 u,v(1,ndim),u3d(1,3),V(1,ndof),kdir(1,ndim),kdir3d(1,3),
    vref(1,ndof),Id_ndim(2,ndim,ndim),tmp1,tmp2,tmp3,tmp4;

  Id_ndim.set(0.);
  for (int jd=1; jd<=ndim; jd++) 
    Id_ndim.setel(1.,jd,jd);

  double p,energy;

  UU.rs().set(U);
  double rho = UU.get(1);
  energy = UU.get(ndim+2)/rho; // total energy
  u.set(UU.is(1,2,ndim+1)).scale(1/rho); UU.is(1);

  double g1=gamma-1;
  static FastMat2 U2,UI2,gU;
  double u2;
  u2 = double(U2.sum_square(u,-1)); 

  double i = energy - 0.5 * u2; // internal energy
  p = g1*rho*i; // pressure
  double cc = sqrt(gamma*p/rho);

  A_jac.rs().set(0.);
  flux.set(0.);
  for (int jdim=1; jdim<=ndim; jdim++) {

    double uk = u.get(jdim);
    // kdir.set(0.).setel(1.,jdim);

    // the corners
    A_jac.rs().ir(1,jdim)
      .setel(uk*(g1*u2-gamma*energy),ndof,1)
      .setel(gamma*uk,ndof,ndof)
      .setel(1.,1,jdim+1);
    
    // vector like submatrices
    A_jac.rs().ir(1,jdim).is(2,2,ndim+1).ir(3,1).axpy(u,-uk)
      .addel(0.5*g1*u2,jdim);
    A_jac.rs().setel(g1,jdim,1+jdim,ndof);
    A_jac.rs().ir(1,jdim).ir(2,ndof).is(3,2,ndim+1).axpy(u,-g1*uk)
      .addel(p/rho+energy,jdim);

    A_jac.rs().ir(1,jdim)
      .is(2,2,ndim+1).is(3,2,ndim+1).set(Id_ndim).scale(uk) // copies uk*identity
      .ir(3,jdim+1).add(u).ir(3) // copies u
      .ir(2,jdim+1).axpy(u,-g1).rs();
    
    flux.rs().ir(2,jdim).prod(A_jac.rs().ir(1,jdim),UU,1,-1,-1);
    flux.rs();
  }

  if ( options & COMP_UPWIND ) {

    // A_grad_U es ndof
    A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);
    Uintri.prod(iJaco,u,1,-1,-1);
    double vel = sqrt(u2);
    UI2.sum_square(Uintri,-1);
    double h_supg = 2.*vel/sqrt(double(UI2));

    double c = sqrt(gamma*p/rho);
    lam_max = c+vel;
    // double lam_max = c;
    tau_a = tau_fac * h_supg/(2.* lam_max);
    double vmax = -1;

    vref.setel(UU.get(1),1);
    vref.is(1,2,ndim+1).set(vel*rho).rs();
    vref.setel(UU.get(ndof),ndof);

    for (int jdof=1; jdof<=ndof; jdof++) {
           double vaux = vref.get(jdof);
	   vaux = vaux*vaux;
	   gU.sum_square(grad_U.rs().ir(2,jdof),-1);
	   vaux = double(gU)/vaux;
	   vmax = (vmax > vaux ? vmax : vaux);
    }
    grad_U.rs();
    vmax = sqrt(vmax);

    // Shock Capturing term. If not used return tau_supg as usual and
    // delta_sc=0.
    tau_delta = 0; delta_sc=0.;
    double rho_i_inv = 1./(rho*i);
#ifdef USE_FASTMAT2_CACHE
    FastMat2::branch();
#endif
    if ( shock_capturing && (vmax > shock_capturing_threshold/h_supg) ) {
#ifdef USE_FASTMAT2_CACHE
      FastMat2::choose(0);
#endif
      // calculo del tensor metrico de Riemann (A0) para transformar de variables 
      // conservativas a variables de entropia (ver Mallet o Shakib)
      V.setel(0.,1); 
      V.is(1,2,ndim+1).set(UU.is(1,2,ndim+1)).scale(rho_i_inv).rs();
      UU.rs();
      double Vndof = -UU.get(1)*rho_i_inv;
      V.setel(Vndof,ndof);

      v.set(V.is(1,2,ndim+1)); V.rs();
      double k1  = 0.5*v.sum_square_all()/Vndof;
      double ZZZ = (k1+1)*Vndof;

      A01.setel(k1*k1+gamma,1,1)
	.setel(ZZZ,1,ndof)
	.setel(ZZZ,ndof,1)
	.setel(Vndof*Vndof,ndof,ndof);
      A01.ir(1,1).is(2,2,ndim+1).set(v).scale(k1).rs();
      A01.ir(2,1).is(1,2,ndim+1).set(v).scale(k1).rs();

      A01.ir(2,ndof).is(1,2,ndim+1).set(v).scale(Vndof).rs();
      A01.ir(1,ndof).is(2,2,ndim+1).set(v).scale(Vndof).rs();

      A01.is(1,2,ndim+1).is(2,2,ndim+1).prod(v,v,1,2);
      for (int jjj=1; jjj<=ndim; jjj++) 
	A01.addel(-Vndof,jjj,jjj);
      A01.rs().scale(1./(rho*i*Vndof));

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

    }
#ifdef USE_FASTMAT2_CACHE
    FastMat2::leave();
#endif

    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
    ret_options |= SCALAR_TAU;
    tau_supg.setel(tau_supg_d,1,1);
  } 

  if (options & COMP_SOURCE) {
    G_source.set(0.);
  }
  return 0;
 
#if 0 // fixme:= falta convertir
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
#endif
}

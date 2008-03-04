//$Id: strat_sw2d.cpp,v 1.8 2007/01/30 19:03:44 rodrigop Exp $
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#include <stdio.h>
#include <string.h>
#include <cassert>
#include <vector>
#include <src/util2.h>
#include <src/utils.h>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "strat_sw2d.h"

// State vector U = [u1 v1 h1 u2 v2 h2]
#define NDOF 6
#undef AJACX
#undef AJACY
#define AJACX(j,k) VEC2(ajacx,((j)-1),((k)-1),NDOF)
#define AJACY(j,k) VEC2(ajacy,((j)-1),((k)-1),NDOF)


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::start_chunk(int &options) {
  options |= SCALAR_TAU;	// tell the advective element routine
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
  
  int ierr;
  int nel,nelprops;
  elemset->elem_params(nel,ndof,nelprops);
  //o Acceleration of gravity.
  EGETOPTDEF_ND(elemset,double,gravity,1.);
  //o rho1
  EGETOPTDEF_ND(elemset,double,rho1,0.);
  //o rho2
  EGETOPTDEF_ND(elemset,double,rho2,0.);
  PETSCFEM_ASSERT0((rho1>0. && rho2>0.),"densities must be positive");
  //o Scale the SUPG upwind term. 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);
  //o Add shock-capturing term.
  EGETOPTDEF_ND(elemset,int,shock_capturing,0);
  // o Add shock-capturing term if relative variation of variables
  // inside the element exceeds this.
  // _T: double[ndim]/double[ndim*ndof]/double[ndim*ndof*ndof] 
  //  _N: advective_jacobians _D: no default  _DOC: 
  //i_tex ../../doc/advdifop.tex advective_jacobians
  //  _END
  //o Threshold value for height.
  EGETOPTDEF_ND(elemset,double,h_min,1e-6);
  //o Threshold value for velocity.
  EGETOPTDEF_ND(elemset,double,vel_min,1e-6);
  assert(ierr==0);
  //o Kinematic viscosity (nu=mu/rho)
  EGETOPTDEF_ND(elemset,double,nu,1.e-5);
  //o Diffusive jacobians factor
  EGETOPTDEF_ND(elemset,double,diff_factor,1.);
  //o Threshold value for $h$ while computing turbulence model.
  EGETOPTDEF_ND(elemset,double,h_min,1e-6);
  
  //o Dimension of the problem. 
  EGETOPTDEF_ND(elemset,int,ndim,0);
  PETSCFEM_ASSERT0(ndim==2,"Only 2D Stratified Shallow Water eqs.");
  A_jac.resize(3,ndim,ndof,ndof);
  D_jac.resize(4,ndim,ndim,ndof,ndof);
  //  flux_mom1.resize(2,ndim,ndim);
  //  flux_mom2.resize(2,ndim,ndim);
  C_jac.resize(2,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  W_N.resize(2,nel,nel);
  tmp1.resize(1,ndof);
  tmp2.resize(2,nel,nel);
  tmp3.resize(2,ndof,ndof);
  vref.resize(1,ndof);
  UU.resize(1,ndof);
  u1m.resize(1,ndim);
  u2m.resize(1,ndim);
  flux_mass1.resize(1,ndim);
  flux_mass2.resize(1,ndim);
  Uintri.resize(1,ndim);
  Uintri1.resize(1,ndim);
  Uintri2.resize(1,ndim);
  A01.resize(2,ndof,ndof);
  bottom_slope.resize(1,ndim);
  grad_U_psi.resize(2,ndim,ndof);
  tmp33.resize(2,ndof,ndim);
  tmp11.resize(4,nel,ndim,ndof,ndof);

  // strat sw 2d must be used with weak form 0 because flux functions
  // are not conservative
  EGETOPTDEF_ND(elemset,int,weak_form,1);
  PETSCFEM_ASSERT0((weak_form==0),"weak_form must be zero for 2D stratified sw eqs.");
#ifdef USE_A_JAC_DUMMY
  //para debug de caso lineal
  A_jac_dummy.resize(3,ndim,ndof,ndof);  
  A_jac_dummy.set(0.).setel(0.,1,1,2).setel(1.,1,2,1);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stratsw2d_ff::stratsw2d_ff(const NewElemset *e) 
  : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
stratsw2d_ff::~stratsw2d_ff() {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::set_state(const FastMat2 &U) {
  UU.rs().set(U);
  h1=UU.get(ndim+1);
  h2=UU.get(2*(ndim+1));
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::enthalpy(FastMat2 &H) {
  H.set(0.);
  double h1_tmp=UU.get(ndim+1);
  double h2_tmp=UU.get(2*(ndim+1));
  H.setel(UU.get(1)*h1_tmp,1);
  H.setel(UU.get(2)*h1_tmp,2);
  H.setel(h1_tmp,3);
  H.setel(UU.get(4)*h2_tmp,4);
  H.setel(UU.get(5)*h2_tmp,5);
  H.setel(h2_tmp,6);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_Cp_N.set(0.);
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
  W_Cp_N.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

void stratsw2d_ff::set_Ufluid(FastMat2 &Uref, FastMat2 &Ufluid) { 
  Ufluid.is(1,1,ndim).set(Uref.rs().is(1,1,ndim));
  Uref.rs();Ufluid.rs();
  Ufluid.is(1,ndim+2,ndim+3).set(Uref.rs().is(1,ndim+2,ndim+3));
  Uref.rs();Ufluid.rs();
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::compute_flux(const FastMat2 &U,
				 const FastMat2 &iJaco, FastMat2 &H,
				 FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
				 FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
				 FastMat2 &tau_supg, double &delta_sc,
				 double &lam_max,FastMat2 &nor, FastMat2 &lambda,
				 FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  static double ajacx[NDOF*NDOF],ajacy[NDOF*NDOF];
  int ierr;
  
  ndof = U.dim(1);
  if ((ndim!=2) || (ndof!=6)) {
    PetscPrintf(PETSC_COMM_WORLD,"Stop stratifief shallow water 2D over 2D domain Only...\n");
    PetscFinalize();
    exit(0);
  }

  static int flag=0;
  static double g=gravity;

  double tau_a, tau_delta, gU, A01v[6*6];
  //  static vector<double> bottom_slope_v;
  
  //  const char *bs;
  //  VOID_IT(bottom_slope_v);
  //   elemset->get_entry("bottom_slope",bs); 
  //   bottom_slope.set(0.);
  //   if (bs!=NULL) {
  //     read_double_array(bottom_slope_v,bs);
  //     assert(bottom_slope_v.size()==(unsigned int)ndim);
  //     bottom_slope.set(&*bottom_slope_v.begin());
  //   }
  
  for (int jj=0; jj<NDOF*NDOF; jj++) {
    ajacx[jj]=0.;
    ajacy[jj]=0.;
  }

  set_state(U);
  flux_mass1.set(UU.is(1,1,ndim));
  UU.rs();
  u1m.set(flux_mass1).scale(1./h1);
  flux_mass2.set(UU.is(1,ndim+2,ndim+3));
  UU.rs();
  u2m.set(flux_mass2).scale(1./h2);

  double uc1 = u1m.sum_square_all();
  double uc2 = u2m.sum_square_all();
  double q1 = sqrt(uc1);
  double q2 = sqrt(uc2);

  double u1,v1,u2,v2;
  u1=u1m.get(1);
  v1=u1m.get(2);
  u2=u2m.get(1);
  v2=u2m.get(2);

  AJACX(1,1)=2*u1;
  AJACX(1,3)=-1.*SQ(u1)+g*h1;
  AJACX(1,6)=g*h1*rho2/rho1;
  AJACX(2,1)=v1;
  AJACX(2,2)=u1;
  AJACX(2,3)=-1.*u1*v1;
  AJACX(3,1)=1.;
  AJACX(4,3)=g*h2;
  AJACX(4,4)=2.*u2;
  AJACX(4,6)=-1.*SQ(u2)+g*h2;
  AJACX(5,4)=v2;
  AJACX(5,5)=u2;
  AJACX(5,6)=-1.*u2*v2;
  AJACX(6,4)=1.;

  A_jac.ir(1,1).set(ajacx);

  AJACY(1,1)=v1;
  AJACY(1,2)=u1;
  AJACY(1,3)=-1.*u1*v1;
  AJACY(2,2)=2.*v1;
  AJACY(2,3)=-1.*SQ(v1)+g*h1;
  AJACY(2,6)=g*h1*rho2/rho1;
  AJACY(3,2)=1.;
  AJACY(4,4)=v2;
  AJACY(4,5)=u2;
  AJACY(4,6)=-1.*u2*v2;
  AJACY(5,3)=g*h2;
  AJACY(5,5)=2*v2;
  AJACY(5,6)=-1.*SQ(v2)+g*h2;
  AJACY(6,5)=1.;

  A_jac.ir(1,2).set(ajacy).rs();

  flux.set(NAN);
  fluxd.set(0.);
  
  //Enthalpy jacobian
  Cp.set(0.);
  Cp.setel(h1,1,1);
  Cp.setel(u1,1,3);
  Cp.setel(h1,2,2);
  Cp.setel(v1,2,3);
  Cp.setel(1.,3,3);
  Cp.setel(h2,4,4);
  Cp.setel(u2,4,6);
  Cp.setel(h2,5,5);
  Cp.setel(v2,5,6);
  Cp.setel(1.,6,6);
  Cp.rs();

  if (options & COMP_UPWIND) {

    D_jac.set(0.);
    C_jac.set(0.);
    // A_grad_U es ndof x 1
    A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);

    Uintri1.prod(iJaco,u1m,1,-1,-1);
    Uintri2.prod(iJaco,u2m,1,-1,-1);
    double h_supg;

    double vel1 = sqrt(uc1);
    double vel2 = sqrt(uc2);


    double vel,h;
    FastMat2::branch();
    if (vel1>=vel2){
      FastMat2::choose(0);
      vel=vel2;
      h=h2;
      Uintri.set(Uintri2);
    } else {
      FastMat2::choose(1);
      vel=vel1;
      h=h1;
      Uintri.set(Uintri1);
    }
    FastMat2::leave();
	
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

    double lam_max = fabs(sqrt((h1+h2)*g)+vel);
    options |= SCALAR_TAU;
    tau_a = SQ(2.*lam_max/h_supg);
    tau_a = tau_fac/sqrt(tau_a);

    double vmax = -1;

    vref.is(1,1,ndim).set(vel*h).rs();
    vref.setel(h,ndim+1).rs();
    vref.is(1,ndim+2,ndim+3).set(vel*h).rs();
    vref.setel(h,ndof).rs();

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
    if (shock_capturing && (vmax > shock_capturing_threshold/h_supg) ) {
      FastMat2::choose(0);
      // calculo del tensor metrico de Riemann (A0) para transformar variables 

      A01v[0]= 1.;
      A01v[1]= 0. ;
      A01v[2]= -u1;

      A01v[3]= 0.;
      A01v[4]= 1.;
      A01v[5]= -v1;

      A01v[6]= -u1;
      A01v[7]= -v1;
      A01v[8]= g*h1+u1*u1+v1*v1;

      // A01.ir(1,1).set(A01v).scale(1./(g*h1)).rs();

      A01v[9]= 1.;
      A01v[10]= 0. ;
      A01v[11]= -u2;

      A01v[12]= 0.;
      A01v[13]= 1.;
      A01v[14]= -v2;

      A01v[15]= -u2;
      A01v[16]= -v2;
      A01v[17]= g*h2+u2*u2+v2*v2;

      // A01.ir(1,2).set(A01v).scale(1./(g*h2)).rs();
      A01.set(A01v).scale(1./(g*(h1+h2))).rs();
      // calculo del delta shock capturing delta_sc
      double vaux_num,vaux_den;
      tmp1.prod(A01,A_grad_U,1,-1,-1);
      tmp22.prod(A_grad_U,tmp1,-1,-1);
      vaux_num = double(tmp22);

      grad_U_psi.prod(iJaco,grad_U,1,-1,-1,2);
      tmp33.prod(A01,grad_U_psi,1,-1,2,-1);
      tmp4.prod(grad_U_psi,tmp33,-1,-2,-2,-1);
      vaux_den = double(tmp4);

      delta_sc = sqrt(vaux_num/vaux_den);
      tau_delta = delta_sc/(lam_max*lam_max);
      
    } // if (shock_capturing ...

    FastMat2::leave();

    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
    options |= SCALAR_TAU;
    tau_supg.setel(tau_supg_d,1,1);

  } 

  if (options & COMP_SOURCE) {
    G_source.set(0.);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(A_jac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,C_jac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,C_jac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void stratsw2d_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
					  FastMat2 &grad_N,double w) {
  grad_N_D_grad_N.set(0.);
  tmp11.prod(D_jac,grad_N,-1,2,3,4,-1,1).scale(w);
  grad_N_D_grad_N.prod(tmp11,grad_N,1,-1,2,4,-1,3);
}

void stratsw2d_ff::get_Cp(FastMat2 &Cp_a) {
  Cp_a.set(Cp);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void stratsw2d_ff::get_Ajac(FastMat2 &Ajac_a) {
  Ajac_a.set(A_jac);
}

void stratsw2d_ff::Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
				FastMat2 &Rie, FastMat2 &drdU,
				FastMat2 &C_U){
#if 0  
//   int ierr;
//   // FIxME:= gravity: Why this is here again? It is
//   // already above...
//   //o Gravity of the problem.
//   EGETOPTDEF_ND(elemset,double,gravity,1.);
//   double tmpd,tmpd1,tmpd2,tmpd3,tt,tt2,pp,
//     ppg1,ppg,h_eps=1.e-10,signudn=0.0,u_eps=1.e-10;
//   // riemann invariants, jacobians and characteristics for sw2d
//   tmpd1=U.get(1); tmpd3=normal.get(1); tmpd2=U.get(ndof);
//   pp=tmpd1*tmpd3;
//   signudn=(tmpd3>0.0 ? 1.0 : -1.0);
//   //  if (fabs(pp)<u_eps) signudn=0.0;
//   tt=2.*sqrt(gravity*tmpd2);
//   Rie.setel(pp+tt,1).setel(pp-tt,ndof);
//   tt=sqrt(gravity*tmpd2);
//   C_U.setel(pp+tt,1).setel(pp-tt,ndof);
//   ppg1=(tmpd2>h_eps ? sqrt(gravity/tmpd2) : 0.);//fix me with area and wl_width
//   ppg=0.;
//   drdU.setel(signudn,1,1).setel(signudn,2,1)
//     .setel(ppg+ppg1,1,2).setel(ppg-ppg1,2,2);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void stratsw2d_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(A_jac,normal,-1,1,2,-1);
}

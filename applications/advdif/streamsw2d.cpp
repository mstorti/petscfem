//$Id: streamsw2d.cpp,v 1.8 2007/01/30 19:03:44 mstorti Exp $
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

#include "streamsw2d.h"

#define NDOF 3
#undef AJACX
#undef AJACY
#define AJACX(j,k) VEC2(ajacx,((j)-1),((k)-1),NDOF)
#define AJACY(j,k) VEC2(ajacy,((j)-1),((k)-1),NDOF)


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::start_chunk(int &options) {
  options |= SCALAR_TAU;	// tell the advective element routine
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
  
  int ierr;
  int nel,nelprops;
  elemset->elem_params(nel,ndof,nelprops);
  //o Acceleration of gravity.
  EGETOPTDEF_ND(elemset,double,gravity,1.);
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
  //o Correcting factor for diffusion in the $k$ transport equation
  EGETOPTDEF_ND(elemset,double,sigma_k,1.);
  //o Correcting factor for diffusion in the $\epsilon$ transport equation
  EGETOPTDEF_ND(elemset,double,sigma_e,1.3);
  //o Coefficient for shallow water turbulent model
  EGETOPTDEF_ND(elemset,double,C_mu,0.09);
  //o Coefficient for shallow water turbulent model
  EGETOPTDEF_ND(elemset,double,C_1,1.44);
  //o Coefficient for shallow water turbulent model
  EGETOPTDEF_ND(elemset,double,C_2,1.92);
  //o Coefficient for shallow water turbulent model
  EGETOPTDEF_ND(elemset,double,D,1.);
  //o Chezy coefficient for bottom friction modelling
  EGETOPTDEF_ND(elemset,double,Chezy,110);
  //o Threshold value for $\epsilon$ (clip below this)
  EGETOPTDEF_ND(elemset,double,eps_min,1e-6);
  //o Threshold value for $k$ while computing turbulence model.
  EGETOPTDEF_ND(elemset,double,ket_min,1e-6);
  //o Threshold value for $h$ while computing turbulence model.
  EGETOPTDEF_ND(elemset,double,h_min,1e-6);
  //o Threshold value for velocity while computing turbulence model.
  EGETOPTDEF_ND(elemset,double,vel_min,1e-6);
  
  //o Dimension of the problem. 
  EGETOPTDEF_ND(elemset,int,ndim,0);
  assert(ndim==2);
  A_jac.resize(3,ndim,ndof,ndof);
  D_jac.resize(4,ndim,ndim,ndof,ndof);
  flux_mom.resize(2,ndim,ndim);
  C_jac.resize(2,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  W_N.resize(2,nel,nel);
  tmp1.resize(1,ndof);
  tmp2.resize(2,nel,nel);
  tmp3.resize(2,ndof,ndof);
  vref.resize(1,ndof);
  u.resize(1,ndim);
  flux_mass.resize(1,ndim);
  Uintri.resize(1,ndim);
  A01.resize(2,ndof,ndof);
  dev_tens.resize(2,2,2);
  tmp5.resize(2,2,2);
  bottom_slope.resize(1,ndim);
  grad_U_psi.resize(2,ndim,ndof);
  tmp33.resize(2,ndof,ndim);
  tmp11.resize(4,nel,ndim,ndof,ndof);
#ifdef USE_A_JAC_DUMMY
  //para debug de caso lineal
  A_jac_dummy.resize(3,ndim,ndof,ndof);  
  A_jac_dummy.set(0.).setel(0.,1,1,2).setel(1.,1,2,1);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
streamsw2d_ff::streamsw2d_ff(const NewElemset *e) 
  : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
streamsw2d_ff::~streamsw2d_ff() {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::set_state(const FastMat2 &U) {
  UU.rs().set(U);
  h=UU.get(ndim+1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::enthalpy(FastMat2 &H) {
  H.set(0.);
  double h_tmp=UU.get(3);
  H.setel(UU.get(1)*h_tmp,1);
  H.setel(UU.get(2)*h_tmp,2);
  H.setel(h_tmp,3);
  //esto esta OK!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_Cp_N.set(0.);
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
  W_Cp_N.rs();
  //esta cuenta esta OK!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
  //esta cuenta esta OK!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::compute_flux(const FastMat2 &U,
				 const FastMat2 &iJaco, FastMat2 &H,
				 FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
				 FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
				 FastMat2 &tau_supg, double &delta_sc,
				 double &lam_max,FastMat2 &nor, FastMat2 &lambda,
				 FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  static double ajacx[NDOF*NDOF],ajacy[NDOF*NDOF];
  int ierr;
  
  if ((ndim!=2) || (ndof!=3)) {
    PetscPrintf(PETSC_COMM_WORLD,"Stop shallow_water 2D over 2D domain Only...\n");
    PetscFinalize();
    exit(0);
  }

  static int flag=0;
  static double g=gravity;

  double tau_a, tau_delta, gU, A01v[9];
  static vector<double> bottom_slope_v;
  ndof = U.dim(1);
  
  const char *bs;
  VOID_IT(bottom_slope_v);
  elemset->get_entry("bottom_slope",bs); 
  bottom_slope.set(0.);
  if (bs!=NULL) {
    read_double_array(bottom_slope_v,bs);
    assert(bottom_slope_v.size()==(unsigned int)ndim);
    bottom_slope.set(&*bottom_slope_v.begin());
  }
    
  for (int jj=0; jj<NDOF*NDOF; jj++) {
    ajacx[jj]=0.;
    ajacy[jj]=0.;
  }

  set_state(U);
  flux_mass.set(UU.is(1,1,ndim));
  UU.rs();
  u.set(flux_mass).scale(1./h);

  double u2 = u.sum_square_all();
  double q = sqrt(u2);

  double ux,uy;
  ux=u.get(1);
  uy=u.get(2);

  AJACX(1,1) = 2*ux;
  AJACX(1,2) = 0.;
  AJACX(1,3) = -ux*ux+g*h;
  AJACX(2,1) = uy;
  AJACX(2,2) = ux;
  AJACX(2,3) = -ux*uy;
  AJACX(3,1) = 1.;
  AJACX(3,2) = 0.;
  AJACX(3,3) = 0.;

  A_jac.ir(1,1).set(ajacx);

  AJACY(1,1) = uy;
  AJACY(1,2) = ux;
  AJACY(1,3) = -ux*uy;
  AJACY(2,1) = 0.;
  AJACY(2,2) = 2*uy;
  AJACY(2,3) = -uy*uy + g*h;
  AJACY(3,1) = 0.;
  AJACY(3,2) = 1.;
  AJACY(3,3) = 0.;

  A_jac.ir(1,2).set(ajacy).rs();

  flux_mom.prod(u,u,1,2).scale(h);

  double h_term = 0.5*g*h*h;
  for (int jdim=1; jdim<=ndim; jdim++) {
    flux_mom.addel(h_term,jdim,jdim);
  }
  flux.rs().is(1,1,ndim).set(flux_mom);
  flux.rs().ir(1,3).set(flux_mass);
  flux.rs();

  //Enthalpy jacobian
  Cp.set(0.);
  Cp.setel(h,1,1);
  Cp.setel(0.,1,2);
  Cp.setel(ux,1,3);
  Cp.setel(0.,2,1);
  Cp.setel(h,2,2);
  Cp.setel(uy,2,3);
  Cp.setel(0.,3,1);
  Cp.setel(0.,3,2);
  Cp.setel(1.,3,3);
  Cp.rs();



  if (options & COMP_UPWIND) {

    D_jac.set(0.);
    double nu_h= nu/h;
    D_jac.setel( 2.*nu_h   ,1,1,1,1);
    D_jac.setel(-2.*ux*nu_h,1,1,1,3);
    
    D_jac.setel(-uy*nu_h   ,1,1,2,3);
    D_jac.setel( nu_h      ,1,1,2,2);
    
    D_jac.setel( nu_h      ,2,2,1,1);
    D_jac.setel(-ux*nu_h   ,2,2,1,3);
    
    D_jac.setel(-2.*uy*nu_h,2,2,2,3);
    D_jac.setel( 2.*nu_h   ,2,2,2,2);
    
    D_jac.setel( nu_h      ,1,2,2,1);
    D_jac.setel(-ux*nu_h   ,1,2,2,3);
    
    D_jac.setel(-uy*nu_h   ,2,1,1,3);
    D_jac.setel( nu_h      ,2,1,1,2);
    
    D_jac.scale(diff_factor);

    grad_U.is(2,1,2);
    dev_tens.set(grad_U);
    tmp5.set(dev_tens).t();
    dev_tens.add(tmp5);

    dev_tens.scale(nu);
    grad_U.rs();
    fluxd.set(0.).is(1,1,2).add(dev_tens).rs();
    fluxd.scale(diff_factor);

    double vel = sqrt(u2);

    // Code C_jac here...

    double tmp61=SQ(Chezy)*(h<h_min ? h_min : h)*(vel<vel_min ? vel_min : vel);
    C_jac.set(0.);

    C_jac.setel(-g*(SQ(ux)+u2)/tmp61,1,1);
    C_jac.setel(-g*ux*uy/tmp61,1,2);
    C_jac.setel(-g*grad_H.get(1,1)+2*g*u2*ux/tmp61,1,3);

    C_jac.setel(-g*(SQ(uy)+u2)/tmp61,2,2);
    C_jac.setel(-g*ux*uy/tmp61,2,1);
    C_jac.setel(-g*grad_H.get(2,1)+2*g*u2*uy/tmp61,2,3);

    C_jac.scale(-1.0);

    // A_grad_U es ndof x 1
    A_grad_U.rs().prod(A_jac.rs(),grad_U,-1,1,-2,-1,-2);

    Uintri.prod(iJaco,u,1,-1,-1);
    double h_supg;

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

    lam_max = fabs(sqrt(h*g)+vel);

    options |= SCALAR_TAU;
    tau_a = SQ(2.*lam_max/h_supg);
    tau_a = tau_fac/sqrt(tau_a);

    double vmax = -1;

    vref.is(1,1,ndim).set(vel*h).rs();
    double pp = UU.get(ndof);
    vref.setel(pp,ndof);

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
      // calculo del tensor metrico de Riemann (A0) para transformar de variables 

      A01v[0]= 1.;
      A01v[1]= 0. ;
      A01v[2]= -ux;

      A01v[3]= 0.;
      A01v[4]= 1.;
      A01v[5]= -uy;

      A01v[6]= -ux;
      A01v[7]= -uy;
      A01v[8]= g*h+ux*ux+uy*uy;

      A01.set(A01v).scale(1./(g*h));

      // calculo del delta shock capturing delta_sc
      double vaux_num,vaux_den;
      tmp1.prod(A01,A_grad_U,1,-1,-1);
      tmp22.prod(A_grad_U,tmp1,-1,-1);
      vaux_num = double(tmp22);

      grad_U_psi.prod(iJaco,grad_U,1,-1,-1,2);
      tmp33.prod(A01,grad_U_psi,1,-1,2,-1);
      tmp4.prod(grad_U_psi,tmp33,-1,-2,-2,-1);
      vaux_den = double(tmp4);

      delta_sc = sqrt(vaux_num / vaux_den);
      tau_delta = delta_sc/(lam_max*lam_max);
      
    } // if (shock_capturing ...

    FastMat2::leave();

    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
    options |= SCALAR_TAU;
    tau_supg.setel(tau_supg_d,1,1);

  } 

  if (options & COMP_SOURCE) {
    grad_H.ir(2,1);
    G_source
      .set(0.)
      .is(1,1,ndim)
      .add(grad_H)
      // .add(bottom_slope)  ??????????
      .scale(-g*h)
      .axpy(u,-g/SQ(Chezy)*q)
      .rs();
    grad_H.rs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(A_jac,grad_N,-1,2,3,-1,1);
  //esta cuenta esta OK!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,C_jac,1,3,2,4);
  //esta cuenta esta OK!!!!
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,C_jac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
  //esto esta OK!!!
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw2d_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
					  FastMat2 &grad_N,double w) {
  grad_N_D_grad_N.set(0.);
  tmp11.prod(D_jac,grad_N,-1,2,3,4,-1,1).scale(w);
  grad_N_D_grad_N.prod(tmp11,grad_N,1,-1,2,4,-1,3);
  //esta cuenta esta OK!!!
}

void streamsw2d_ff::Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
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
//   /* 
//      if ((area<1.e-6) || (wl_width<1.e-7)) {
//      tt=0.0;
//      ppg=0.0;
//      ppg1=0.0;
//      tt2=0.0;
//      } else {
//      tt=2.*sqrt(gravity*area/wl_width);
//      ppg=(tmpd1/area)*(wl_width-1.);
//      ppg1=sqrt(gravity/(wl_width*area));
//      tt2=sqrt(gravity*area/wl_width);
//      }
//   */
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void streamsw2d_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(A_jac,normal,-1,1,2,-1);
  //esta cuenta esta OK!!!
}

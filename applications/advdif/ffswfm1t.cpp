//__INSERT_LICENSE__
//$Id: ffswfm1t.cpp,v 1.1 2003/04/01 12:16:14 mstorti Exp $

#include <stdio.h>
#include <string.h>
#include <vector>
#include <cassert>

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/utils.h>

#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define NDOF 5
#define AJACX(j,k) VEC2(ajacx,((j)-1),((k)-1),NDOF)
#define AJACY(j,k) VEC2(ajacy,((j)-1),((k)-1),NDOF)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "swfm2t_ff_t::operator()"
int swfm2t_ff_t::operator()(ADVDIFFF_ARGS) {

  static double ajacx[NDOF*NDOF],ajacy[NDOF*NDOF];
  int ierr;

  if (ndim!=2) {
    PFEMERRQ("Stop with \"mojito\" drinking!! 2D Only...\n");
  }

  static int flag=0,ndof;
  static double g, gravity, tau_fac,shock_capturing_threshold,
    sigma_k,sigma_e,C_mu,C_1,C_2,D,Chezy,C_P_e,eps_min,ket_min,
    P_h,P_e,P_k,h_min,vel_min;
  double tau_a, tau_delta, gU, A01v[9];
  static int shock_capturing;
  static vector<double> bottom_slope_v;
  static FastMat2 bottom_slope(1,ndim);

  // Load properties only once.

  FastMat2::branch();
  if (start_chunk) {
    FastMat2::choose(0);
    start_chunk = 0;

//    if (start_chunk) {
//      start_chunk = 1;

    ndof = U.dim(1);

    C_jac.resize(2,ndof,ndof);

    //o Acceleration of gravity
    EGETOPTDEF_ND(elemset,double,gravity,1.);
    g=gravity;
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
    const char *bs;
    VOID_IT(bottom_slope_v);
    elemset->get_entry("bottom_slope",bs); 
    bottom_slope.set(0.);
    if (bs!=NULL) {
      read_double_array(bottom_slope_v,bs);
      assert(bottom_slope_v.size()==ndim);
      bottom_slope.set(bottom_slope_v.begin());
    }
    
    
    for (int jj=0; jj<NDOF*NDOF; jj++) {
      ajacx[jj]=0.;
      ajacy[jj]=0.;
    }
    AJACX(3,1) = 1.;
    AJACY(3,2) = 1.;

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
    EGETOPTDEF_ND(elemset,double,h_min,1e-6);

    C_P_e = C_2*sqrt(C_mu)*pow(g,1.25)/sqrt(D)/pow(Chezy,2.5);
  }
  FastMat2::leave();

  static FMatrix u(ndim),flux_mass(ndim),UU,flux_mom,
    Uintri,vref(ndof),A01(ndof,ndof),tmp1,tmp2,tmp3,tmp4,grad_U_psi,
    dev_tens(2,2),tmp5,tmp61,tmp62,tmp63,tmp64,tmp65,tmp66;
  UU.rs().set(U);

  double h = UU.get(ndim+1);
  double HH = H.get(1);
  flux_mass.set(UU.is(1,1,ndim));
  UU.rs();
  u.set(flux_mass).scale(1./h);
  double ket = UU.get(4)/h;
  double eps = UU.get(5)/h;

  // Turbulent kinematic viscosity
  double nu_t = C_mu * SQ(ket) / (eps<eps_min ? eps_min : eps);

  double u2 = u.sum_square_all();
  double q = sqrt(u2);

  double ux,uy;
  ux=u.get(1);
  uy=u.get(2);

  AJACX(1,1) = 2*ux;
  AJACX(1,3) = -ux*ux+g*h;
  AJACX(2,1) = uy;
  AJACX(2,2) = ux;
  AJACX(2,3) = -ux*uy;
  AJACX(4,4) = ux;
  AJACX(5,5) = ux;

  A_jac.ir(1,1).set(ajacx);

  AJACY(1,1) = uy;
  AJACY(1,2) = ux;
  AJACY(1,3) = -ux*uy;
  AJACY(2,2) = 2*uy;
  AJACY(2,3) = -uy*uy + g*h;
  AJACY(4,4) = uy;
  AJACY(5,5) = uy;

  A_jac.ir(1,2).set(ajacy).rs();

  flux_mom.prod(u,u,1,2).scale(h);

  // double h_term = 0.5*g*(h*h-HH*HH); !! ERROR
  double h_term = 0.5*g*h*h-(2./3.)*ket*h; 
  for (int jdim=1; jdim<=ndim; jdim++) {
    flux_mom.addel(h_term,jdim,jdim);
  }
  flux.rs().is(1,1,ndim).set(flux_mom);
  flux.rs().ir(1,3).set(flux_mass);
  flux.rs().ir(1,4).set(flux_mass).scale(ket);
  flux.rs().ir(1,5).set(flux_mass).scale(eps);
  flux.rs();


  if (options & COMP_UPWIND) {

    // Code D_jac here...
    D_jac.set(0.);
    double nu_t_h= nu_t / h;
    D_jac.setel( 2.*nu_t_h   ,1,1,1,1);
    D_jac.setel(-2.*ux*nu_t_h,1,1,1,3);

    D_jac.setel(-uy*nu_t_h   ,1,1,2,3);
    D_jac.setel( nu_t_h      ,1,1,2,2);

    D_jac.setel( nu_t/sigma_k,1,1,4,4);
    D_jac.setel( nu_t/sigma_e,1,1,5,5);

    D_jac.setel( nu_t_h      ,2,2,1,1);
    D_jac.setel(-ux*nu_t_h   ,2,2,1,3);

    D_jac.setel(-2.*uy*nu_t_h,2,2,2,3);
    D_jac.setel( 2.*nu_t_h   ,2,2,2,2);

    D_jac.setel( nu_t/sigma_k,2,2,4,4);
    D_jac.setel( nu_t/sigma_e,2,2,5,5);

    D_jac.setel( nu_t_h      ,1,2,2,1);
    D_jac.setel(-ux*nu_t_h   ,1,2,2,3);

    D_jac.setel(-uy*nu_t_h   ,2,1,1,3);
    D_jac.setel( nu_t_h      ,2,1,1,2);


    // Turbulent viscous stresses
    grad_U.is(2,1,2);
    dev_tens.set(grad_U);
    tmp5.set(dev_tens).t();
    dev_tens.add(tmp5);

    // Production of kinetic turbulence
    double tmp62 = 0.5*dev_tens.sum_square_all();

    P_h = tmp62*nu_t/h;

    dev_tens.scale(nu_t);
    grad_U.rs();
    fluxd.set(0.).is(1,1,2).add(dev_tens).rs();

    P_k = g*CB(q)/SQ(Chezy);
    P_e = C_P_e*pw4(q)/h;

    // Turbulent diffusion of k and e
    grad_U.ir(2,4);
    fluxd.ir(1,4).axpy(grad_U,nu_t/sigma_k);
    grad_U.ir(2,5);
    fluxd.ir(1,5).axpy(grad_U,nu_t/sigma_e);
    grad_U.rs(); fluxd.rs();

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

    double tmp63=C_2*sqrt(C_mu)*pow(g,1.25)/sqrt(D)/pow(Chezy,2.5);
    double tmp64=C_1*C_mu*tmp62/SQ((h<h_min?h_min:h));
    double tmp65=SQ((h<h_min?h_min:h))*(eps<eps_min?eps_min:eps);
    double tmp66=SQ((h<h_min?h_min:h))*SQ((eps<eps_min?eps_min:eps));

    C_jac.setel(3*g*u2*ux/tmp61,4,1);
    C_jac.setel(3*g*u2*uy/tmp61,4,2);
    C_jac.setel(-2*C_mu*SQ(ket)/tmp65*tmp62-3*g*SQ(u2)/tmp61,4,3);
    C_jac.setel( 2*C_mu*ket/tmp65*tmp62,4,4);
    C_jac.setel(-C_mu*SQ(ket)/tmp66*tmp62-1.,4,5);

    C_jac.setel(4*tmp63*ux*u2/SQ((h<h_min?h_min:h)),5,1);
    C_jac.setel(4*tmp63*uy*u2/SQ((h<h_min?h_min:h)),5,2);
    C_jac.setel(-2*tmp64*ket-5*tmp63*SQ(u2)/SQ((h<h_min?h_min:h)),5,3);
    C_jac.setel(tmp64+C_2*SQ(eps)/SQ((ket<ket_min?ket_min:ket)),5,4);
    C_jac.setel(-2*C_2*eps/(ket<ket_min?ket_min:ket),5,5);

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

    ret_options |= SCALAR_TAU;
    tau_a = SQ(2.*lam_max/h_supg)
	  +18.*SQ(nu_t/SQ(h_supg));
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
    if ( shock_capturing && (vmax > shock_capturing_threshold /h_supg) ) {
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
      tmp2.prod(A_grad_U,tmp1,-1,-1);
      vaux_num = double(tmp2);

      grad_U_psi.prod(iJaco,grad_U,1,-1,-1,2);
      tmp3.prod(A01,grad_U_psi,1,-1,2,-1);
      tmp4.prod(grad_U_psi,tmp3,-1,-2,-2,-1);
      vaux_den = double(tmp4);

      delta_sc = sqrt(vaux_num / vaux_den);
      tau_delta = delta_sc/(lam_max*lam_max);
      
    } // if (shock_capturing ...

    FastMat2::leave();

    double tau_supg_d = ((tau_a-tau_delta)>0 ? (tau_a-tau_delta) : 0);
    ret_options |= SCALAR_TAU;
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
    G_source.setel(P_h+P_k-h*eps,4);
    G_source.setel((C_1*P_h-C_2*eps*h)*eps/(ket<ket_min ? ket_min : ket)+P_e,5);
  }
}

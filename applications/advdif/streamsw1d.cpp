//$Id: streamsw1d.cpp,v 1.3.2.2 2004/02/25 14:27:34 mstorti Exp $
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#include <stdio.h>
#include <string.h>
#include <cassert>
#include <src/util2.h>
#include <src/utils.h>
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "streamsw1d.h"

//NDOF=2:(hu h)
#define NDOF 2
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#undef AJAC
#define AJAC(j,k) VEC2(ajac,((j)-1),((k)-1),NDOF)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::start_chunk(int &options) {
  options |= SCALAR_TAU;	// tell the advective element routine
  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
  if (!channel) channel = ChannelShape::factory(elemset);
  channel->init();
  
  if (!friction_law) 
    friction_law = FrictionLaw::factory(elemset);
  friction_law->init();
  
  int ierr;
  int nel,nelprops;
  ndimel=1;
  elemset->elem_params(nel,ndof,nelprops);
  // Acceleration of gravity
  EGETOPTDEF_ND(elemset,double,gravity,1.);
  // Scale the SUPG upwind term. 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);
  // Threshold value for height 
  EGETOPTDEF_ND(elemset,double,h_min,1e-6);
  // Threshold value for velocity
  EGETOPTDEF_ND(elemset,double,vel_min,1e-6);
  // Scale friction term
  EGETOPTDEF_ND(elemset,double,cfric,1.);
  assert(ierr==0);

  EGETOPTDEF_ND(elemset,int,ndim,0);
  assert(ndim==2);

  A_jac.resize(3,ndimel,ndof,ndof);
  //  flux_mass.resize(2,ndimel,ndimel);
  flux_mom.resize(2,ndimel,ndimel);
  C_jac.resize(2,ndof,ndof);
  Sf_jac.resize(2,ndof,ndof);
  Cp.resize(2,ndof,ndof);
  W_N.resize(2,nel,nel);
  tmp2.resize(2,nel,nel);
  tmp3.resize(2,ndof,ndof);
  // vref.resize(1,ndof);
#ifdef USE_A_JAC_DUMMY
  //para debug de caso lineal
  A_jac_dummy.resize(3,ndimel,ndof,ndof);  
  A_jac_dummy.set(0.).setel(0.,1,1,2).setel(1.,1,2,1);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::element_hook(ElementIterator &element) {
  channel->element_hook(element);
  friction_law->element_hook(element);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
streamsw1d_ff::streamsw1d_ff(const NewElemset *e) 
  : AdvDifFFWEnth(e), channel(NULL), 
  friction_law(NULL) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
streamsw1d_ff::~streamsw1d_ff() { delete friction_law; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::set_state(const FastMat2 &U) {
  UU.rs().set(U);
  double h =UU.get(ndof);
  channel->geometry(h,area,wl_width,perimeter);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::enthalpy(FastMat2 &H) {
  H.set(0.);
  H.setel(area*UU.get(1),1);
  H.setel(area,2);
  H.scale(tmp_mask);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_Cp_N.set(0.);
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
  W_Cp_N.rs();
  /*
    W_Cp_N.ir(2,1).ir(4,1);
    W_Cp_N.prod(W,N,1,2).scale(weight);
    W_Cp_N.ir(2,2).ir(4,2);
    W_Cp_N.prod(W,N,1,2).scale(weight);
    W_Cp_N.rs();
  */
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/* void streamsw1d_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
   //  P_Cp.set(P_supg).scale(wl_width);
   } */

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  delta_sc = 1.;		// For smoothing 
  adv_mask = 1.;
  tmp_mask = 1.;

  // In quasilinear form: dU/dt + A*dU/dx = G
  // A:= the flux jacobian dF/dU 
  // We use primitive variables in present context: [u h]
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  static double ajac[NDOF*NDOF];
  if ((ndim!=2) || (dim()!=1) || (ndof!=2)) {
    PetscPrintf(PETSC_COMM_WORLD,"Stop shallow_water 1D over 2D domain Only...\n");
    PetscFinalize();
    exit(0);
  }
  double tau_a, tau_delta, gU, g=gravity;
  double Sf;//c_Sf_jac_1,c_Sf_jac_2;
  // Load properties only once.
  set_state(U,grad_U);
  for (int jj=0; jj<ndof*ndof; jj++) {
    ajac[jj]=0.;
  }

  double h=UU.get(ndof);//es la comp 2 de U
  double u=UU.get(1);//
  channel->geometry(h,area,wl_width,perimeter);

  friction_law->flow_Sf(area,perimeter,u,Sf,Sf_jac);

  double h_tmp=(h<h_min ? h_min : h);//para terminos de friccion y gravedad

  flux_mass.set(UU.is(1,1,1).scale(area));
  UU.rs();
  
  double u2=u*u;
  double q=sqrt(u2);
  
  double ux;
  ux=u;
  
  AJAC(1,1)=2*ux*area;
  AJAC(1,2)=(ux*ux+0.5*g*h)*wl_width+0.5*g*area; 
  AJAC(2,1)=area;
  AJAC(2,2)=ux*wl_width;
  
  A_jac.ir(1,1).set(ajac).rs();
  A_jac.scale(adv_mask);

  flux_mom.setel(area*ux*ux,1,1);
  // le agrego el termino con h
  double h_term=0.5*g*area*h;
  flux_mom.addel(h_term,1,1);
  
  flux.rs().is(1,1,1).set(flux_mom);
  flux.rs().ir(1,2).set(flux_mass);
  flux.rs();
  flux.scale(adv_mask);

#ifdef USE_A_JAC_DUMMY
  flux.prod(A_jac_dummy,UU,2,1,-1,-1);
  A_jac.set(A_jac_dummy);
#endif

  // Si no es difusivo hay que ponerlo a 0!!!!
  fluxd.set(0.);

  //Enthalpy jacobian
  Cp.set(0.);
  Cp.setel(area,1,1);
  //  Cp.setel(h*wl_width,1,2);
  Cp.setel(ux*wl_width,1,2);
  Cp.setel(0.,2,1);
  Cp.setel(wl_width,2,2);
  Cp.rs();
  Cp.scale(tmp_mask);
  
  if (options & COMP_UPWIND) {
    
    double vel=sqrt(u2);
    
    // Code C_jac here... that's for Newton loop of the right hand side of s-w
    C_jac.set(0.);
    C_jac.setel(cfric*(-g*(Sf_jac.ir(1,1).get(1))),1,1);
    //    Sf_jac.rs();
    double gUtmp1=-0.5*g*grad_U.ir(1,1).get(2); grad_U.rs();
    double gUtmp2=cfric*(g*Sf_jac.get(ndof));
    double gUtmp3=g*grad_H.get(1,1);
    C_jac.setel(wl_width*(gUtmp1+gUtmp2+gUtmp3),1,2);
    C_jac.scale(-1.0);
    //C_jac.scale(.0); //se lo pongo solo por canal rectangular!!!!!!!!!!! fix this!!

    Sf_jac.rs();
    // A_grad_U es ndof x 1
    A_grad_U.prod(A_jac,grad_U,-1,1,-2,-1,-2);
    
    Uintri.set(iJaco).scale(u);//iJaco es la inv del jac de la transf de cordenadas
    double h_supg;
    
    FastMat2::branch();
    if (vel>1e-10) {
      FastMat2::choose(0);
      h_supg = 2.*vel/sqrt(Uintri.sum_square_all());
    } else {
      FastMat2::choose(1);
      h_supg = 2./sqrt(iJaco.sum_square_all());
    }
    FastMat2::leave();
    
    lam_max = fabs(sqrt(h*g)+vel); //max eigenvalue 
#ifdef USE_A_JAC_DUMMY
    lam_max = 1.;
#endif
    
    tau_a = SQ(2.*lam_max/h_supg);
    tau_a = tau_fac/sqrt(tau_a);
    tau_supg.setel(tau_a,1,1);
  } 

  double pq=grad_U.ir(1,1).get(2); grad_U.rs(); 
  double ppq=grad_H.get(1,1); grad_H.rs();
  if (options & COMP_SOURCE) {
    G_source
      .set(0.)
      .is(1,1,ndimel)
      .add(-g*area*ppq)
      .add(-0.5*g*pq*(h*wl_width-area))
      .add(cfric*g*area*Sf) //le saco el signo menos, supuestamente es asi
      .rs();
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(A_jac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,C_jac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,C_jac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void streamsw1d_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & grad_N,double w) {
  grad_N_D_grad_N.set(0.);
}

void streamsw1d_ff::Riemann_Inv(const FastMat2 &U, const FastMat2 &normaln,
				FastMat2 &Rie, FastMat2 &drdU,
				FastMat2 &C_U){
  int ierr;
  EGETOPTDEF_ND(elemset,double,gravity,1.);
  double tmpd,tmpd1,tmpd2,tmpd3,tt,tt2,pp,
    ppg1,ppg,h_eps=1.e-10,signudn=0.0,u_eps=1.e-10;
  // riemann invariants, jacobians and characteristics for sw1d
  tmpd1=U.get(1); tmpd3=normaln.get(1); tmpd2=U.get(ndof);
  pp=tmpd1*tmpd3;
  signudn=(tmpd3>0.0 ? 1.0 : -1.0);
  //  if (fabs(pp)<u_eps) signudn=0.0;
  tt=2.*sqrt(gravity*tmpd2);
  Rie.setel(pp+tt,1).setel(pp-tt,ndof);
  tt=sqrt(gravity*tmpd2);
  C_U.setel(pp+tt,1).setel(pp-tt,ndof);
  ppg1=(tmpd2>h_eps ? sqrt(gravity/tmpd2) : 0.);//fix me with area and wl_width
  ppg=0.;
  drdU.setel(signudn,1,1).setel(signudn,2,1)
    .setel(ppg+ppg1,1,2).setel(ppg-ppg1,2,2);
  /* 
     if ((area<1.e-6) || (wl_width<1.e-7)) {
     tt=0.0;
     ppg=0.0;
     ppg1=0.0;
     tt2=0.0;
     } else {
     tt=2.*sqrt(gravity*area/wl_width);
     ppg=(tmpd1/area)*(wl_width-1.);
     ppg1=sqrt(gravity/(wl_width*area));
     tt2=sqrt(gravity*area/wl_width);
     }
  */
}

//__INSERT_LICENSE__
//$Id: gaschem.cpp,v 1.7 2003/11/13 02:49:39 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "gaschem.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::start_chunk(int &ret_options) {
  int ierr;

  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
  elemset->elem_params(nel,ndof,nelprops);
  //o Dimension of the problem
  EGETOPTDEF_ND(elemset,int,ndim,0);
  //o Gas constant
  EGETOPTDEF_ND(elemset,double,Rgas,8.314);
  //o Gas temperature
  EGETOPTDEF_ND(elemset,double,Tgas,293.0);
  //o Direction of gravity
  EGETOPTDEF_ND(elemset,int,g_dir,ndim);
  //o Turbulent viscosity
  EGETOPTDEF_ND(elemset,double,nu_t,0.);
  assert(nu_t>0.);
  //o Turbulent Schmidt number
  EGETOPTDEF_ND(elemset,double,Sc,0.83);
  //o Scales bubble/liquid film exchange
  EGETOPTDEF_ND(elemset,double,hm_fac,1.0);

  //o Cutoff value for Nb
  EGETOPTDEF_ND(elemset,double,Nb_ctff,0.0);
  //o Cutoff value for CO
  EGETOPTDEF_ND(elemset,double,CO_ctff,0.0);
  //o Cutoff value for CN
  EGETOPTDEF_ND(elemset,double,CN_ctff,0.0);
  //o Cutoff value for CdO
  EGETOPTDEF_ND(elemset,double,CdO_ctff,0.0);
  //o Cutoff value for CdN
  EGETOPTDEF_ND(elemset,double,CdN_ctff,0.0);
  //o Scale Nb eq. with this
  EGETOPTDEF_ND(elemset,double,Nb_scale,1.0);

  // Values for KO and KN from [Buscaglia et.al. 2002]
  //are in mol/m3/bar

  //o Henry constant for O2 #[mol/m3/Pa]#
  EGETOPTDEF_ND(elemset,double,KO,1.3516e-5);
  //o Henry constant for N2 #[mol/m3/Pa]#
  EGETOPTDEF_ND(elemset,double,KN,0.6788e-5);

  U.resize(1,ndof);
  Cp.resize(2,ndof,ndof);
  Cp.eye();
  Ajac.resize(3,ndim,ndof,ndof);
  Ajac.set(0.);
  u_liq.resize(1,ndim);
  u_gas.resize(1,ndim);

  Djac.resize(4,ndim,ndof,ndim,ndof).set(0.);
  Cjac.resize(2,ndof,ndof).set(0.);
  Uintri.resize(1,ndim);
  tmp0.resize(1,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
gaschem_ff::gaschem_ff(NewElemset *e) : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
gaschem_ff::~gaschem_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::set_state(const FastMat2 &UU) {
  U.set(UU);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::enthalpy(FastMat2 &H) {
  H.set(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
#if 0
void gaschem_ff::compute_tau(int ijob) {
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
// compute the profile for each specific subproblem according to jobinfo value
void gaschem_ff::set_profile(FastMat2 &seed) {
  seed.eye();
  // seed.set(1);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau

  H.is(1,1,ndim);
  u_liq.set(H);
  H.rs();

  double diff;
#define GC_VAR(name,indx) double name = ctff(U.get(indx),diff,name##_ctff);
  GC_VAR(Nb,1);
  GC_VAR(CO,2);
  GC_VAR(CN,3);
  GC_VAR(CdO,4);
  GC_VAR(CdN,5);
  Nb = Nb*Nb_scale;
  
  double pgas = H.get(ndim+1);
  double vb = (CO+CN)*Rgas*Tgas/(pgas*Nb);
  double rb = pow(vb*3.0/(4.0*M_PI),(1.0/3.0));
  double vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
		  rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
  u_gas.set(u_liq).addel(vslip,g_dir);
  
  // Nb, CO and CN are advected with u_gas
  for (int j=1; j<=3; j++) Ajac.ir(2,j).ir(3,j).set(u_gas);
  // CdO and CdN are advected with u_liq
  for (int j=4; j<=5; j++) Ajac.ir(2,j).ir(3,j).set(u_liq);
  Ajac.rs();

  A_grad_U.prod(Ajac,grad_U,-1,1,-2,-1,-2);

  for (int j=1; j<=3; j++) flux.ir(1,j).set(u_gas).scale(U.get(j));
  for (int j=4; j<=5; j++) flux.ir(1,j).set(u_liq).scale(U.get(j));
  flux.rs();

  double Dg = nu_t/Sc;
  for (int j=1; j<=5; j++) {
    grad_U.ir(2,j);
    fluxd.ir(1,j).set(grad_U).scale(Dg);
  }
  grad_U.rs();
  fluxd.rs();

  for (int j=1; j<=5; j++) Djac.ir(2,j).ir(4,j).eye(Dg);
  Djac.rs();

  if (options & COMP_UPWIND) {
    advdf_e = new_adv_dif_elemset;
#define pi M_PI
    double Volume = advdf_e->volume();
    int axi = advdf_e->axi;
    double h_pspg;

    if (ndim==2 | (ndim==3 && axi>0)) {
      h_pspg = sqrt(4.*Volume/pi);
    } else if (ndim==3) {
      h_pspg = cbrt(6*Volume/pi);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Only dimensions 2 and 3 allowed for this element.\n");
    }

    Uintri.prod(iJaco,u_gas,1,-1,-1);

    // this is approx. 2*U/h
    double Uh = sqrt(Uintri.sum_square_all());
    double tau, tau_fac=1.0;
    double vel = sqrt(u_gas.sum_square_all());
    lam_max = vel;

    FastMat2::branch();
    if (vel*vel > 1e-5*Uh*Dg) { // remove singularity when v=0
      FastMat2::choose(0);
      double Pe  = vel*vel/(Uh*Dg);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = tau_fac/Uh*magic; // intrinsic time
    } else {
      FastMat2::choose(1);
      double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
      tau = tau_fac*h*h/(12.*Dg);
    }
    FastMat2::leave();

    tau_supg.d(1,2).is(1,1,3).set(tau);
    
    Uintri.prod(iJaco,u_liq,1,-1,-1);

    // this is approx. 2*U/h
    Uh = sqrt(Uintri.sum_square_all());
    vel = sqrt(u_liq.sum_square_all());
    if (vel>lam_max) lam_max = vel;

    FastMat2::branch();
    if (vel*vel > 1e-5*Uh*Dg) { // remove singularity when v=0
      FastMat2::choose(0);
      double Pe  = vel*vel/(Uh*Dg);	// Peclet number
      // magic function
      double magic = (fabs(Pe)>1.e-4 ? 1./tanh(Pe)-1./Pe : Pe/3.); 
      tau = tau_fac/Uh*magic; // intrinsic time
    } else {
      FastMat2::choose(1);
      double h = 2./sqrt(tmp0.sum_square(iJaco,1,-1).max_all());
      tau = tau_fac*h*h/(12.*Dg);
    }
    FastMat2::leave();
    tau_supg.is(1,4,5).set(tau).rs();

  }

  if (options & COMP_SOURCE) {
    // Source terms ================
    // Film coefficient
    double rb_crit = 6.67e-4;
    double hm = hm_fac*(rb < rb_crit ? (rb/rb_crit)*4e-4 : 4e-4);
    
    double coef = 3.0*(CO+CN)*Rgas*Tgas*hm/(pgas*rb);
    double xO = CO/(CO+CN);
    double xN = 1.0-xO;
    
    double SO = coef*(CdO-KO*pgas*xO);
    double SN = coef*(CdN-KN*pgas*xN);
    
    G_source.setel(0.,1);
    G_source.setel(SO,2);
    G_source.setel(SN,3);
    G_source.setel(-SO,4);
    G_source.setel(-SN,5);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(Ajac,normal,-1,1,2,-1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(Ajac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				     FastMat2 &dshapex,double w) {
  tmp1.prod(Djac,dshapex,-1,2,3,4,-1,1).scale(w);
  grad_N_D_grad_N.prod(tmp1,dshapex,1,2,-1,4,-1,3);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,Cjac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void gaschem_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,Cjac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}

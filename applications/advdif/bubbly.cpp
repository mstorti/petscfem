//__INSERT_LICENSE__
//$Id: bubbly.cpp,v 1.9 2002/02/28 22:29:21 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::start_chunk(int &ret_options) {
  int ierr;
  FastMat2 tmp5;

  elemset->elem_params(nel,ndof,nelprops);  
  EGETOPTDEF_ND(elemset,int,ndim,0);
  EGETOPTDEF_ND(elemset,double,rho_l,0);
  EGETOPTDEF_ND(elemset,double,rho_g,0);
  EGETOPTDEF_ND(elemset,double,visco_l,0);
  EGETOPTDEF_ND(elemset,double,visco_g,0);

  // Turbulence parameters
  //o C_mu
  EGETOPTDEF_ND(elemset,double,C_mu,0.09);
  //o C_1
  EGETOPTDEF_ND(elemset,double,C_1,1.44);
  //o C_2
  EGETOPTDEF_ND(elemset,double,C_2,1.92);
  //o sigma_k 
  EGETOPTDEF_ND(elemset,double,sigma_k,1.);
  //o sigma_e 
  EGETOPTDEF_ND(elemset,double,sigma_e,1.3);
  //o sigma_e 
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  //o Adjust the stability parameters, taking into account
  // the time step. If the \verb+steady+ option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  // \verb+temporal_stability_factor+ is set to 0.
  EGETOPTDEF_ND(elemset,double,temporal_stability_factor,1.);

  //o _T: double[ndim] _N: G_body _D: null vector 
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.resize(1,ndim);
  G_body.set(0.);
  ierr = elemset->get_double("G_body",
			     *G_body.storage_begin(),1,ndim);

  assert(ndim>0);
  assert(ndof==4+2*ndim);
  assert(rho_l>0.);
  assert(visco_l>0.);
  assert(rho_g>0.);
  assert(visco_g>0.);
  vl_indx = 3;
  vl_indxe = 3+ndim-1;
  vg_indx = vl_indx+ndim;
  vg_indxe = vg_indx+ndim-1;
  k_indx = vg_indx+ndim;
  e_indx = k_indx+1;
  v_l.resize(1,ndim);
  v_g.resize(1,ndim);
  v_mix.resize(1,ndim);
  Cp.resize(2,ndof,ndof);
  Cpc.resize(2,ndof,ndof);
  Ajac.resize(3,ndim,ndof,ndof);
  Ajacc.resize(3,ndim,ndof,ndof);
  Id.resize(2,ndim,ndim).eye();
  Amoml.resize(2,ndim,ndim);
  Amomg.resize(2,ndim,ndim);
  Y.resize(3,ndim,ndim,ndim);

  Djac.resize(4,ndim,ndof,ndim,ndof);
  tmp1.resize(4,nel,ndof,ndim,ndof);
  Cjac.resize(2,ndof,ndof);
  Cjacc.resize(2,ndof,ndof);
  tmp2.resize(2,nel,nel);
  tmp3.resize(2,ndof,ndof);

  grad_v_l.resize(2,ndim,ndim);
  strain_rate_l.resize(2,ndim,ndim);
  strain_rate_g.resize(2,ndim,ndim);
  grad_v_g.resize(2,ndim,ndim);
  grad_k.resize(1,ndim);
  grad_e.resize(1,ndim);
  
  IdId.resize(4,ndim,ndim,ndim,ndim);
  IdId.prod(Id,Id,1,3,2,4);
  tmp5.prod(Id,Id,2,3,1,4);
  IdId.add(tmp5);

  uintri.resize(1,ndim);
  svec.resize(1,ndim);
  tmp9.resize(1,nel);

  Djacc.resize(4,ndim,ndof,ndim,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bubbly_ff::bubbly_ff(const NewAdvDif *e) 
  : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
bubbly_ff::~bubbly_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::set_state(const FastMat2 &UU) {
  U.set(UU);

  // Scalar variables
  p = U.get(1);
  alpha_l = U.get(2);
  alpha_g = 1.-alpha_l;
  arho_l = alpha_l*rho_l;
  arho_g = alpha_g*rho_g;

  k = U.get(k_indx);
  eps = U.get(e_indx);
  // Velocities
  U.is(1,vl_indx,vl_indxe);
  v_l.set(U);
  U.rs().is(1,vg_indx,vg_indxe);
  v_g.set(U);
  U.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::set_state(const FastMat2 &U,const FastMat2 &grad_U) {
  set_state(U);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::enthalpy(FastMat2 &H) {
  // Scalar values
  H.setel(arho_l,1);
  H.setel(arho_g,2);
  H.setel(arho_l*k,k_indx);
  H.setel(arho_l*eps,e_indx);
  // Vector values
  H.is(1,vl_indx,vl_indxe).set(v_l).scale(arho_l);
  H.rs();
  H.is(1,vg_indx,vg_indxe).set(v_g).scale(arho_g);
  H.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff
::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  double strain_rate_scalar;

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Enthalpy Jacobian Cp
  // First column
  Cp.set(0.);
  Cp.setel(rho_l,1,1);
  Cp.setel(-rho_g,2,1);
  Cp.is(1,vl_indx,vl_indxe).ir(2,1).set(v_l).scale(rho_l);
  Cp.is(1).is(1,vg_indx,vg_indxe).set(v_g).scale(-rho_g);
  Cp.rs();
  Cp.setel(rho_l*k,k_indx,1);
  Cp.setel(rho_l*eps,e_indx,1);
  // V_l column block
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe)
    .eye(arho_l);
  // V_g column block
  Cp.rs().is(1,vg_indx,vg_indxe).is(2,vg_indx,vg_indxe)
    .eye(arho_g);
  // Turbulence
  Cp.rs();
  Cp.setel(arho_l,k_indx,k_indx);
  Cp.setel(arho_l,e_indx,e_indx);

  Cpc.set(Cp).is(2,2).is(2,1).is(2,3,ndof);
  Cp.set(Cpc);
  Cpc.rs();
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Advective fluxes
  flux.set(0.);
  flux.ir(1,1).set(v_l).scale(arho_l);
  flux.ir(1,2).set(v_g).scale(arho_g);

  flux.rs();
  Amoml.prod(v_l,v_l,1,2).scale(rho_l)
    .axpy(Id,p);
  flux.is(1,vl_indx,vl_indxe).axpy(Amoml,alpha_l);

  flux.rs();
  Amomg.prod(v_g,v_g,1,2).scale(rho_g)
    .axpy(Id,p);
  flux.is(1,vg_indx,vg_indxe).axpy(Amomg,alpha_g);

  // TUrb
  flux.rs();
  flux.ir(1,k_indx).set(v_l).scale(arho_l*k);
  flux.ir(1,e_indx).set(v_l).scale(arho_l*eps);
  flux.rs();
  // flux.set(0.);  //debug:= 
  // Cambiamos signo de la ec. de cont. debug:=
  flux.is(1,1).scale(-1.).rs(); 
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // Adjective Jacobians
  // First column
  Ajac.set(0.);
  Ajac.ir(2,1).ir(3,1).set(v_l).scale(rho_l);
  Ajac.ir(2,2).set(v_g).scale(-rho_g);
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).set(Amoml);
  Ajac.rs().is(2,vg_indx,vg_indxe).ir(3,1).axpy(Amomg,-1.);
  Ajac.rs().ir(2,k_indx).ir(3,1).set(v_l).scale(rho_l*k);
  Ajac.rs().ir(2,e_indx).ir(3,1).set(v_l).scale(rho_l*eps);
  
  // Second column
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,2).axpy(Id,alpha_l);
  Ajac.rs().is(2,vg_indx,vg_indxe).ir(3,2).axpy(Id,alpha_g);
  
  Ajac.rs().ir(2,1).is(3,vl_indx,vl_indxe).axpy(Id,arho_l);
  Ajac.rs().ir(2,2).is(3,vg_indx,vg_indxe).axpy(Id,arho_g);
  
  /// fixme:= Verify this!!!
  Y.prod(v_l,Id,2,1,3).scale(arho_l);
  Ajac.rs().is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe)
    .add(Y);
  Y.prod(v_l,Id,1,2,3).scale(arho_l);
  Ajac.add(Y);

  Y.rs().prod(v_g,Id,1,2,3).scale(arho_g);
  Ajac.rs().is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe)
    .add(Y);
  Y.exc(2,3);
  Ajac.add(Y);
  Y.rs();
  
  Ajac.rs().ir(2,k_indx).is(3,vl_indx,vl_indxe)
    .axpy(Id,arho_l*k);
  Ajac.rs().ir(2,e_indx).is(3,vl_indx,vl_indxe)
    .axpy(Id,arho_l*eps);

  Ajac.rs().ir(2,k_indx).ir(3,k_indx).set(v_l).scale(arho_l);
  Ajac.rs().ir(2,e_indx).ir(3,e_indx).set(v_l).scale(arho_l);
  Ajac.rs();

  Ajacc.set(Ajac).is(3,2).is(3,1).is(3,3,ndof);
  Ajac.set(Ajacc);
  Ajacc.rs();
  // Ajac.set(0.); //debug:= 
  // Cambiamos signo de la ec. de cont. debug:=
  Ajac.is(2,1).scale(-1.).rs(); 

  A_grad_U.prod(Ajac,grad_U,-1,1,-2,-1,-2);

  // Strain rate for the liquid
  // fixme:= not clear which indices are... But resulting matrix is
  // symmetric so that it should not matter
  grad_U.is(2,vl_indx,vl_indxe);
  grad_v_l.set(grad_U);
  grad_U.rs();
  strain_rate_l.set(grad_v_l);
  grad_v_l.t();
  strain_rate_l.add(grad_v_l).scale(0.5);
  grad_v_l.rs();

  // Turbulent viscosity
  // fixme:= la tasa de deformacion se toma la del liquido??
  strain_rate_scalar = strain_rate_l.sum_square_all();
  // fixme:= esto no compila 
  visco_t = C_mu*rho_l * square(k)/eps;
  P_k = 2*visco_t*strain_rate_scalar;

  visco_l_eff = visco_l + visco_t;
  visco_g_eff = alpha_g * visco_g + alpha_l * visco_l_eff;
  fluxd.set(0.);

  // Strain rate for the gas
  grad_U.is(2,vg_indx,vg_indxe);
  grad_v_g.set(grad_U);
  grad_U.rs();
  strain_rate_g.set(grad_v_g);
  grad_v_g.t();
  strain_rate_g.add(grad_v_g).scale(0.5);
  grad_v_g.rs();

  fluxd.is(1,vl_indx,vl_indxe).set(strain_rate_l).scale(2.*alpha_l*visco_l_eff);
  fluxd.rs().is(1,vg_indx,vg_indxe).set(strain_rate_g).scale(2.*visco_g_eff);
  fluxd.rs();

  // Turbulent diffusion for k-e
  grad_U.ir(2,k_indx);
  grad_k.set(grad_U);
  fluxd.ir(1,k_indx).set(grad_k).scale(alpha_l*visco_l_eff/sigma_k);
  
  grad_U.ir(2,e_indx);
  grad_e.set(grad_U);
  fluxd.ir(1,e_indx).set(grad_e).scale(alpha_l*visco_l_eff/sigma_e);
  fluxd.rs();
  grad_U.rs();
  
  // Diffusive jacobians
  Djac.set(0.);
  Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,alpha_l*visco_l_eff);
  Djac.rs().is(2,vg_indx,vg_indxe).is(4,vg_indx,vg_indxe).axpy(IdId,visco_g_eff);
  Djac.rs().ir(2,k_indx).ir(4,k_indx).set(Id).scale(alpha_l*visco_l_eff/sigma_k);
  Djac.rs().ir(2,e_indx).ir(4,e_indx).set(Id).scale(alpha_l*visco_l_eff/sigma_e);
  Djac.rs();

  Djacc.set(Djac).is(4,2).is(4,1).is(4,3,ndof);
  Djac.set(Djacc);
  Djacc.rs();
  // Cambiamos signo de la ec. de cont. debug:=
  // Djac.is(2,1).scale(-1.).rs(); // lo sacamos porque no tiene difusion!!!
  
  // Reactive terms
  Cjac.set(0.);
  // fixme:= verificar en las dos sig. lineas que es la primera columna
  Cjac.is(1,vl_indx,vl_indxe).ir(2,1).set(G_body).scale(rho_l);
  Cjac.rs().is(1,vg_indx,vg_indxe).ir(2,1).set(G_body).scale(-rho_g);
  Cjac.rs();

  Cjac.setel(4.*rho_l*C_mu*strain_rate_scalar*k/eps,k_indx,k_indx);
  double Cke = (-2.*rho_l*C_mu*strain_rate_scalar*square(k)+rho_l*square(eps))/square(eps);
  Cjac.setel(-Cke,k_indx,e_indx);
  double Cek = (-2.*rho_l*C_1*strain_rate_scalar*square(k)+C_2*rho_l*square(eps))/square(k);
  Cjac.setel(-Cek,e_indx,k_indx);
  Cjac.setel(-2*eps*C_2*rho_l/k,e_indx,e_indx);
  // fixme:= No sabemos esto bien con que signo va... 
  Cjac.scale(-1.);

  Cjacc.set(Cjac).is(2,2).is(2,1).is(2,3,ndof);
  Cjac.set(Cjacc);
  Cjacc.rs();
  
  if (options & COMP_UPWIND) {
    const NewAdvDif *advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);
#define pi M_PI
    double Volume = advdf_e->volume();
    double h_pspg,Delta;
    if (ndim==2) {
      h_pspg = sqrt(4.*Volume/pi);
      Delta = sqrt(Volume);
    } else if (ndim==3) {
      // h_pspg = pow(6*Volume/pi,1./3.);
      // El pow() da segmentation violation cuando corro con -O !!
      h_pspg = cbrt(6*Volume/pi);
      Delta = cbrt(Volume);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Only dimensions 2 and 3 allowed for this element.\n");
    }

    // fixme:= Deberia estar pesado por las densidades ???
    // v_mix.set(v_l).scale(alpha_l).axpy(v_g,alpha_g);
    // fixme:= Por ahora ponemos v_mix = v_l_o
    // (velocidad del liquido en el paso anterior) 
    
    FastMat2 &Uo = (FastMat2 &) advdf_e->Uold();
    Uo.is(1,vl_indx,vl_indxe);
    v_mix.set(Uo);
    Uo.rs();
    uintri.prod(iJaco,v_mix,1,-1,-1);
    double Uh = uintri.sum_square_all();
    Uh = sqrt(Uh)/2;
    double velmod = sqrt(v_mix.sum_square_all());

    double tol=1.0e-16;
    double h_supg=0;
    const FastMat2 &grad_N = *advdf_e->grad_N();
    FastMat2::branch();
    if(velmod>tol) {
      FastMat2::choose(0);
      // svec:= a streamline oriented unit vector
      svec.set(v_mix).scale(1./velmod);
      h_supg = tmp9.prod(grad_N,svec,-1,1,-1).sum_abs_all();
      h_supg = (h_supg < tol ? tol : h_supg);
      h_supg = 2./h_supg;
    } else {
      h_supg = h_pspg;
    }
    FastMat2::leave();

    // printf("h_supg: %f, h_pspg: %f\n",h_supg,h_pspg);
    // h_supg = h_pspg;
    double Peclet = velmod * h_supg / (2. * visco_l_eff);
    double rec_Dt = advdf_e->rec_Dt();
    double tsf = temporal_stability_factor;
    if (rec_Dt==0.) tsf = 0.;
    double tau_supg_a =  tsf * square(2.*rec_Dt)+square(2.*velmod/h_supg)
      +9.*square(4.*visco_l_eff/square(h_supg));
    tau_supg_a = 1./sqrt(tau_supg_a);

    double pspg_advection_factor=1.,pspg_factor=1.;
    double tau_pspg = tsf * square(2.*rec_Dt)
      +square(pspg_advection_factor*2.*velmod/h_pspg)
      +9.*square(4.*visco_l_eff/square(h_pspg));
    tau_pspg = pspg_factor/sqrt(tau_pspg);

    double fz = (Peclet < 3. ? Peclet/3. : 1.);
    double delta_supg = 0.5*h_supg*velmod*fz;
	
    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }
    tau_supg.eye(tau_supg_a).setel(tau_pspg,1,1).setel(tau_pspg,2,2);
    // tau_supg.eye(0.025); // debug:= 
  }
  if (options & COMP_SOURCE) {
    G_source.set(0.);
    G_source.is(1,vl_indx,vl_indxe).set(G_body).scale(arho_l);
    G_source.rs().is(1,vg_indx,vg_indxe).set(G_body).scale(arho_g);
    G_source.rs()
      .setel(P_k-rho_l*eps,k_indx)
      .setel(eps/k*(C_1*P_k-C_2*rho_l*eps),e_indx);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  assert(0); // fixme:= not defined yet (for use with
	     // absorbing b.c's and bcconv)
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N) {
  A_grad_N.prod(Ajac,grad_N,-1,2,3,-1,1);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				     FastMat2 &dshapex,double w) {
  tmp1.prod(Djac,dshapex,-1,2,3,4,-1,1).scale(w);
  grad_N_D_grad_N.prod(tmp1,dshapex,1,2,-1,4,-1,3);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w) {
  tmp2.prod(N,N,1,2).scale(w);
  N_N_C.prod(tmp2,Cjac,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void bubbly_ff::comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
			   FastMat2 &N,double w) {
  tmp3.prod(P_supg,Cjac,1,-1,-1,2).scale(w);
  N_P_C.prod(tmp3,N,1,3,2);
}


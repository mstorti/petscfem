//__INSERT_LICENSE__
//$Id: bubbly.cpp,v 1.17 2003/09/15 01:17:59 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly.h"

extern const char * jobinfo_fields;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::start_chunk(int &ret_options) {
  int ierr;
  FastMat2 tmp5;

  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
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

  // bubble diameter
  EGETOPTDEF_ND(elemset,double,d_bubble,0);
  EGETOPTDEF_ND(elemset,int,comp_interphase_terms,0);

  //o _T: double[ndim] _N: G_body  _D: <null-vector>
  //   _DOC: Vector of gravity acceleration (must be constant). _END
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
  vel_supg.resize(1,ndim);
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
  grad_v_g.resize(2,ndim,ndim);
  strain_rate_g.resize(2,ndim,ndim);
  grad_k.resize(1,ndim);
  grad_e.resize(1,ndim);
  grad_alpha_g.resize(1,ndim);
  grad_p.resize(1,ndim);

  IdId.resize(4,ndim,ndim,ndim,ndim);
  IdId.prod(Id,Id,1,3,2,4);
  tmp5.prod(Id,Id,2,3,1,4);
  IdId.add(tmp5);

  tmp6.resize(3,nel,ndim,ndim);
  tmp7.resize(2,ndim,nel);

  uintri.resize(1,ndim);
  svec.resize(1,ndim);
  tmp9.resize(1,nel);
  W_N.resize(2,nel,nel);

  Djacc.resize(4,ndim,ndof,ndim,ndof);

  v_g_l.resize(1,ndim);
  Phi_1.resize(2,ndim,ndim);
  Phi_2.resize(2,ndim,ndim);

  id_liq = 1, id_gas = -1;

  v_l_old.resize(1,ndim);
  v_g_old.resize(1,ndim);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::element_hook(ElementIterator &element) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
bubbly_ff::bubbly_ff(NewElemset *e) : AdvDifFFWEnth(e) {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
bubbly_ff::~bubbly_ff() { }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::set_state(const FastMat2 &UU) {
  U.set(UU);

  // PetscScalar variables
  p = U.get(1);

  // modif Sept 2002
  //alpha_l = U.get(2);
  alpha_g = U.get(2);
  //alpha_g = 1.-alpha_l;
  alpha_l = 1.-alpha_g;

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

  int comp_liq_prof  = !strcmp(jobinfo_fields,"liq");
  int comp_gas_prof  = !strcmp(jobinfo_fields,"gas");
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  H.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  H.setel(arho_l,1);
  H.is(1,vl_indx,vl_indxe).set(v_l).scale(arho_l);
  H.rs();
  }
  if (comp_kep_prof) {
  // k-eps turbulence
  H.setel(arho_l*k,k_indx);
  H.setel(arho_l*eps,e_indx);
  }
  if (comp_gas_prof) {
  // Gas phase
  H.setel(arho_g,2);
  H.is(1,vg_indx,vg_indxe).set(v_g).scale(arho_g);
  H.rs();
  }

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff
::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
//void bubbly_ff::compute_tau(const FastMat2 &vel_supg, double visco_supg, const FastMat2 &grad_N,
//double h_pspg,double &tau_supg_a,double &tau_pspg, double &delta_supg, int ijob,
//	 double &velmod, double &h_supg) {

void bubbly_ff::compute_tau(int ijob) {

  //    static FastMat2 svec,tmp9;

    const FastMat2 &grad_N = *advdf_e->grad_N();

    double tol=1.0e-16;
    h_supg=0;
    //    const FastMat2 &grad_N = *advdf_e->grad_N();
    velmod = sqrt(vel_supg.sum_square_all());
    FastMat2::branch();
    if(velmod>tol) {
      FastMat2::choose(0);
      // svec:= a streamline oriented unit vector
      svec.set(vel_supg).scale(1./velmod);
      h_supg = tmp9.prod(grad_N,svec,-1,1,-1).sum_abs_all();
      h_supg = (h_supg < tol ? tol : h_supg);
      h_supg = 2./h_supg;
    } else {
      h_supg = h_pspg;
    }
    FastMat2::leave();

    double Peclet = velmod * h_supg / (2. * visco_supg);
    //    double rec_Dt = advdf_e->rec_Dt();
    // double tsf = temporal_stability_factor;
    // if (rec_Dt==0.) tsf = 0.;
    //    tau_supg_a =  tsf * square(2.*rec_Dt)+square(2.*velmod/h_supg)
    //      +9.*square(4.*visco_supg/rho_supg/square(h_supg));
    tau_supg_a =  square(2.*velmod/h_supg)
      +9.*square(4.*visco_supg/square(h_supg));
    tau_supg_a = 1./sqrt(tau_supg_a);
    double pspg_factor=1.;
    tau_pspg = pspg_factor*tau_supg_a;

    if (ijob==0) {
    double fz = (Peclet < 3. ? Peclet/3. : 1.);
    delta_supg = 0.5*h_supg*velmod*fz;
    } else {
    delta_supg = square(2.*velmod/h_supg);
    delta_supg = (delta_supg < tol ? tol : delta_supg);
    delta_supg = 1./sqrt(delta_supg);
    }

	}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
// compute the profile for each specific subproblem according to jobinfo value
void bubbly_ff::set_profile(FastMat2 &seed) {

  int comp_liq_prof  = !strcmp(jobinfo_fields,"liq");
  int comp_gas_prof  = !strcmp(jobinfo_fields,"gas");
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  //  Matrix seed;
  //  seed= Matrix(ndof,ndof);
  //  seed=0;
  seed.set(0.);

  if (comp_liq_prof) {
    int ip[] = {1,3,4,5};    
    for (int j=0; j<=ndim; j++) {
      for (int k=0; k<=ndim; k++) {
	//        seed(ip[j],ip[k])=1;
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  } else if (comp_gas_prof) {
    int ip[] = {2,2+ndim+1,2+ndim+2,2+ndim+3};
    for (int j=0; j<=ndim; j++) {
      for (int k=0; k<=ndim; k++) {
	//        seed(ip[j],ip[k])=1;
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  } else if (comp_kep_prof) {
    int ip[] = {2+2*ndim+1,2+2*ndim+2};
    for (int j=0; j<=1; j++) {
      for (int k=0; k<=1; k++) {
	//        seed(ip[j],ip[k])=1;
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  }
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_flux(const FastMat2 &U,
	       const FastMat2 &iJaco, FastMat2 &H,
	       FastMat2 &grad_H, FastMat2 &flux, FastMat2 &fluxd,
	       FastMat2 &A_grad_U,FastMat2 &grad_U, FastMat2 &G_source,
	       FastMat2 &tau_supg, double &delta_sc,
	       double &lam_max,FastMat2 &nor, FastMat2 &lambda,
	       FastMat2 &Vr, FastMat2 &Vr_inv,int options) {

  double strain_rate_scalar;

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau

  int comp_liq_prof  = !strcmp(jobinfo_fields,"liq");
  int comp_gas_prof  = !strcmp(jobinfo_fields,"gas");
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Enthalpy Jacobian Cp

  Cp.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Cp.setel(-rho_l,1,1);
  Cp.is(1,vl_indx,vl_indxe).ir(2,1).set(v_l).scale(-rho_l).rs();
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe)
    .eye(arho_l).rs();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  Cp.setel(-rho_l*k,k_indx,1);
  Cp.setel(-rho_l*eps,e_indx,1);
  Cp.setel(arho_l,k_indx,k_indx);
  Cp.setel(arho_l,e_indx,e_indx);
  }

  if (comp_gas_prof) {
  // Gas phase
  Cp.setel(rho_g,2,1);
  Cp.is(1,vg_indx,vg_indxe).ir(2,1).set(v_g).scale(rho_g).rs();
  Cp.rs().is(1,vg_indx,vg_indxe).is(2,vg_indx,vg_indxe)
    .eye(arho_g).rs();
  }

  Cpc.set(Cp).is(2,2).is(2,1).is(2,3,ndof);
  Cp.set(Cpc);
  Cpc.rs();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Advective fluxes
  flux.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  flux.ir(1,1).set(v_l).scale(arho_l).rs();
  Amoml.prod(v_l,v_l,1,2).scale(rho_l)
    .axpy(Id,p);
  flux.is(1,vl_indx,vl_indxe).axpy(Amoml,alpha_l).rs();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  flux.ir(1,k_indx).set(v_l).scale(arho_l*k);
  flux.ir(1,e_indx).set(v_l).scale(arho_l*eps).rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  flux.ir(1,2).set(v_g).scale(arho_g).rs();
  Amomg.prod(v_g,v_g,1,2).scale(rho_g);
  flux.is(1,vg_indx,vg_indxe).axpy(Amomg,alpha_g).rs();
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Adjective Jacobians
  Ajac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Ajac.ir(2,1).ir(3,1).set(v_l).scale(-rho_l);
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).axpy(Amoml,-1);
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,2).axpy(Id,alpha_l);
  Ajac.rs().ir(2,1).is(3,vl_indx,vl_indxe).axpy(Id,arho_l);

  Y.rs().prod(v_l,Id,2,1,3).scale(arho_l);
  Ajac.rs().is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe)
    .add(Y);
  Y.prod(v_l,Id,1,2,3).scale(arho_l).rs();
  Ajac.add(Y);
  Y.rs();
  Ajac.rs();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  Ajac.rs().ir(2,k_indx).ir(3,1).set(v_l).scale(-rho_l*k);
  Ajac.rs().ir(2,e_indx).ir(3,1).set(v_l).scale(-rho_l*eps);
  Ajac.rs().ir(2,k_indx).is(3,vl_indx,vl_indxe)
    .axpy(Id,arho_l*k);
  Ajac.rs().ir(2,e_indx).is(3,vl_indx,vl_indxe)
    .axpy(Id,arho_l*eps);
  Ajac.rs().ir(2,k_indx).ir(3,k_indx).set(v_l).scale(arho_l);
  Ajac.rs().ir(2,e_indx).ir(3,e_indx).set(v_l).scale(arho_l);
  Ajac.rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  Ajac.rs().ir(2,2).ir(3,1).set(v_g).scale(rho_g);
  Ajac.rs().is(2,vg_indx,vg_indxe).ir(3,1).set(Amomg);
  Ajac.rs().ir(2,2).is(3,vg_indx,vg_indxe).axpy(Id,arho_g);

  Y.prod(v_g,Id,1,2,3).scale(arho_g);
  Ajac.rs().is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe)
    .add(Y);
  Y.exc(2,3);
  Ajac.add(Y);
  Y.rs();
  Ajac.rs();

  }

  Ajacc.set(Ajac).is(3,2).is(3,1).is(3,3,ndof);
  Ajac.set(Ajacc);
  Ajacc.rs();

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

  // limito la viscosidad por debajo para evitar valores negativos
  double tol=1e-5;
  double alpha_g_ctf = (alpha_g < tol ? tol : alpha_g);
  alpha_g_ctf = (alpha_g_ctf > (1.0-tol) ? (1.0-tol) : alpha_g_ctf);
  double alpha_l_ctf = 1.0-alpha_g_ctf;

  visco_l_eff = alpha_l_ctf * (visco_l + visco_t);
  visco_g_eff = alpha_g_ctf * visco_g +  visco_l_eff;

  // Strain rate for the gas
  grad_U.is(2,vg_indx,vg_indxe);
  grad_v_g.set(grad_U);
  grad_U.rs();
  strain_rate_g.set(grad_v_g);
  grad_v_g.t();
  strain_rate_g.add(grad_v_g).scale(0.5);
  grad_v_g.rs();

  fluxd.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  fluxd.is(1,vl_indx,vl_indxe).set(strain_rate_l).scale(2.*visco_l_eff).rs();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  grad_U.ir(2,k_indx);
  grad_k.set(grad_U);
  fluxd.ir(1,k_indx).set(grad_k).scale(visco_l_eff/sigma_k);
  grad_U.ir(2,e_indx);
  grad_e.set(grad_U);
  fluxd.ir(1,e_indx).set(grad_e).scale(visco_l_eff/sigma_e);
  fluxd.rs();
  grad_U.rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  fluxd.is(1,vg_indx,vg_indxe).set(strain_rate_g).scale(2.*visco_g_eff).rs();
  }


  Djac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,visco_l_eff).rs();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  Djac.rs().ir(2,k_indx).ir(4,k_indx).set(Id).scale(visco_l_eff/sigma_k).rs();
  Djac.rs().ir(2,e_indx).ir(4,e_indx).set(Id).scale(visco_l_eff/sigma_e).rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  Djac.is(2,vg_indx,vg_indxe).is(4,vg_indx,vg_indxe).axpy(IdId,visco_g_eff).rs();
  }

  Djacc.set(Djac).is(4,2).is(4,1).is(4,3,ndof);
  Djac.set(Djacc);
  Djacc.rs();

  // Reactive terms
  //   [1] bouyancy terms
  //   [2] Phase interaction forces
  Cjac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Cjac.is(1,vl_indx,vl_indxe).ir(2,1).set(G_body).scale(-rho_l).rs();
  }
  if (comp_gas_prof) {
  // Gas phase
  Cjac.is(1,vg_indx,vg_indxe).ir(2,1).set(G_body).scale(rho_g).rs();
  }


  if (comp_interphase_terms==1) {
  assert(d_bubble>0);
  C1_drag = 3./4.*visco_l/d_bubble/d_bubble;
  v_g_l.set(v_g).rest(v_l);
  v_slip = v_g_l.sum_square_all();
  v_slip = sqrt(v_slip);
  Rey_bubble =rho_l*d_bubble*v_slip/visco_l;

    FastMat2::branch();
    if (Rey_bubble>0){
      FastMat2::choose(0);

	  // user defined C_drag function (a pata por ahora)
      // Lai & Salcudean
	  C_drag_ff= 24.0/Rey_bubble;
	  dCDdRe_ff= -24.0/Rey_bubble/Rey_bubble;

	  dRedU = rho_l*d_bubble/visco_l;

	  tmp1_drag = C1_drag*Rey_bubble*C_drag_ff;

          tmp2_drag = tmp1_drag*(1.-2.*alpha_g);

          if (comp_liq_prof) {
          // Liquid phase
	  Cjac.is(1,vl_indx,vl_indxe).ir(2,1).axpy(v_g_l,id_liq*tmp2_drag).rs();
          }
          if (comp_gas_prof) {
          // Gas phase
	  Cjac.is(1,vg_indx,vg_indxe).ir(2,1).axpy(v_g_l,id_gas*tmp2_drag).rs();
          }

          tmp4_drag = C_drag_ff+Rey_bubble*dCDdRe_ff;

          Phi_1.prod(v_g_l,v_g_l,1,2).scale(-tmp4_drag*dRedU/v_slip);
          Phi_2.set(Id).scale(Rey_bubble*C_drag_ff);
          Phi_1.rest(Phi_2);
          Phi_1.scale(C1_drag*alpha_l*alpha_g);

          if (comp_liq_prof) {
          // Liquid phase
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).axpy(Phi_1,id_liq*id_liq).rs();
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vg_indx,vg_indxe).axpy(Phi_1,id_liq*id_gas).rs();
          }
          if (comp_gas_prof) {
          // Gas phase
	  Cjac.is(1,vg_indx,vg_indxe).is(2,vg_indx,vg_indxe).axpy(Phi_1,id_gas*id_gas).rs();
	  Cjac.is(1,vg_indx,vg_indxe).is(2,vl_indx,vl_indxe).axpy(Phi_1,id_gas*id_liq).rs();
	  }
    }
    FastMat2::leave();
  }

  if (comp_kep_prof) {
  // k-eps turbulence
  Cjac.setel(4.*rho_l*C_mu*strain_rate_scalar*k/eps,k_indx,k_indx);
  double Cke = (-2.*rho_l*C_mu*strain_rate_scalar*square(k)+rho_l*square(eps))/square(eps);
  Cjac.setel(-Cke,k_indx,e_indx);
  double Cek = (-2.*rho_l*C_1*strain_rate_scalar*square(k)+C_2*rho_l*square(eps))/square(k);
  Cjac.setel(-Cek,e_indx,k_indx);
  Cjac.setel(-2*eps*C_2*rho_l/k,e_indx,e_indx);
  }

  Cjac.scale(-1.);
  Cjacc.set(Cjac).is(2,2).is(2,1).is(2,3,ndof);
  Cjac.set(Cjacc);
  Cjacc.rs();

  if (options & COMP_UPWIND) {
    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);
#define pi M_PI
    double Volume = advdf_e->volume();
    int axi = advdf_e->axi;

    if (ndim==2 | (ndim==3 && axi>0)) {
      h_pspg = sqrt(4.*Volume/pi);
    } else if (ndim==3) {
      h_pspg = cbrt(6*Volume/pi);
    } else {
      PetscPrintf(PETSC_COMM_WORLD,
		  "Only dimensions 2 and 3 allowed for this element.\n");
    }

    FastMat2 &Uo = (FastMat2 &) advdf_e->Uold();

    Uo.is(1,vl_indx,vl_indxe);
    v_l_old.set(Uo);
    Uo.rs();

    Uo.is(1,vg_indx,vg_indxe);
    v_g_old.set(Uo);
    Uo.rs();

    const FastMat2 &grad_N = *advdf_e->grad_N();

    //    double tau_supg_a,tau_pspg,delta_supg,visco_supg,velmod, h_supg;
    // Fase liquida
    vel_supg.set(v_l_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }

    visco_supg = visco_l_eff/rho_l;
    int ijob =0;
    //    compute_tau(vel_supg, visco_supg, grad_N, h_pspg, tau_supg_a, tau_pspg, delta_supg, ijob,
    // velmod, h_supg);
    compute_tau(ijob);

    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }
    tau_supg.eye(tau_supg_a).setel(delta_supg,1,1);

    // Fase gas
    vel_supg.set(v_g_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }
    visco_supg = visco_g_eff/rho_g;
    ijob=1;
    //    compute_tau(vel_supg, visco_supg, grad_N, h_pspg, tau_supg_a, tau_pspg, delta_supg, ijob,
    //velmod, h_supg);
    compute_tau(ijob);

    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }

    for (int k=vg_indx; k<vg_indxe+1; k++) {
    tau_supg.setel(tau_supg_a,k,k);
    }
    tau_supg.setel(delta_supg,2,2);

    tau_supg_c.set(tau_supg);

    grad_U.ir(2,2);
    grad_alpha_g.set(grad_U);
    grad_U.rs();

    if (comp_gas_prof) {
    // Gas phase

    if (1) {
    // Fixme:> agregamos flujo difusivo numerico a la ecuacion de masa de gas
    //         para ver si es el motivo de la falta de convergencia del esquema
    //    fluxd.ir(1,2).set(grad_alpha_l).scale(-0.5*h_supg*velmod*rho_g);
    //    fluxd.ir(1,2).set(grad_alpha_l).scale(-0.5*h_supg*rho_g);
    //    double div_vg = double(tmp10.ctr(grad_v_g,-1,-1));
    grad_v_g.d(1,2);
    tmp10.sum(grad_v_g,-1);
    grad_v_g.rs();

    double div_vg = double(tmp10);

    //    double Dnum =  (div_vg > 0 ? 0 : -div_vg);
    //    double Dnum = fabs(div_vg);
    double Dnum = velmod*h_supg/2.;
    // double Dnum = velmod;
    //    double Dnum_factor = 0.25;
    double Dnum_factor = 0.10;
    //Dnum = -Dnum*rho_g*Dnum_factor*h_supg*h_supg;
    Dnum = -Dnum*rho_g*Dnum_factor;
    // debug
    fluxd.ir(1,2).set(grad_alpha_g).scale(-Dnum);
    //fluxd.ir(1,2).set(grad_alpha_l).scale(Dnum);
    // end debug
    fluxd.rs();
    //    Djac.ir(2,2).ir(4,2).set(Id).scale(-0.5*h_supg*velmod*rho_g);
    //    Djac.ir(2,2).ir(4,2).set(Id).scale(-0.5*h_supg*rho_g);
    // debug
    Djac.ir(2,2).ir(4,2).set(Id).scale(-Dnum);
    //Djac.ir(2,2).ir(4,2).set(Id).scale(Dnum);
    // end debug
    Djac.rs();
    }
    }
  }

  if (options & COMP_SOURCE) {
    G_source.set(0.);
    // Bouyancy forces
    if (comp_liq_prof) {
    // Liquid phase
    G_source.is(1,vl_indx,vl_indxe).set(G_body).scale(arho_l).rs();
    }
    if (comp_gas_prof) {
    // Gas phase
    G_source.is(1,vg_indx,vg_indxe).set(G_body).scale(arho_g).rs();
    }

  if (comp_interphase_terms==1) {
    FastMat2::branch();
    // Phase interaction forces
    if (Rey_bubble>0){
      FastMat2::choose(0);
	// user defined C_drag function (a pata por ahora)
	tmp3_drag = tmp1_drag*alpha_l*alpha_g;
        if (comp_liq_prof) {
        // Liquid phase
	G_source.is(1,vl_indx,vl_indxe).axpy(v_g_l,id_liq*tmp3_drag).rs();
	}
        if (comp_gas_prof) {
        // Gas phase
	G_source.is(1,vg_indx,vg_indxe).axpy(v_g_l,id_gas*tmp3_drag).rs();
	}
    }
    FastMat2::leave();
  }

    if (comp_kep_prof) {
    // k-eps turbulence
    G_source.setel(P_k-rho_l*eps,k_indx)
            .setel(eps/k*(C_1*P_k-C_2*rho_l*eps),e_indx);
    }

    // DEBUG fixme after try this change
    if (0) {

    // - p * grad (alfa_g) to do the scheme pressure invariant
    grad_U.ir(2,2);
    grad_alpha_g.set(grad_U);
    grad_U.rs();
    G_source.is(1,vl_indx,vl_indxe).axpy(grad_alpha_g,-p).rs();
    G_source.is(1,vg_indx,vg_indxe).axpy(grad_alpha_g,p).rs();
    } else {

    // pressure term to do the system pressure invariant
    if (comp_liq_prof) {
    // Liquid phase
    grad_U.ir(2,2);
    grad_alpha_g.set(grad_U);
    grad_U.rs();
    G_source.is(1,vl_indx,vl_indxe).axpy(grad_alpha_g,-p).rs();
    }
    if (comp_gas_prof) {
    // Gas phase
    grad_U.ir(2,1);
    grad_p.set(grad_U);
    grad_U.rs();
    G_source.is(1,vg_indx,vg_indxe).axpy(grad_p,-alpha_g).rs();
    }
    }
    // end of DEBUG
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
  A_jac_n.prod(Ajac,normal,-1,1,2,-1);
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


#ifdef USE_COMP_P_SUPG
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::comp_P_supg(FastMat2 &P_supg) {

  int comp_liq_prof  = !strcmp(jobinfo_fields,"liq");
  int comp_gas_prof  = !strcmp(jobinfo_fields,"gas");
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  double rho_m,tau;

    const FastMat2 &grad_N = *new_adv_dif_elemset->grad_N();
    // esto es equivalente a lo viejo
    if (0) {
    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    P_supg.prod(Ao_grad_N,tau_supg_c,1,2,-1,-1,3);
    }
    // lo nuevo a ver si anda ????
    // P_supg es de (nel,ndof,ndof)
    else {
    P_supg.set(0.);
    if (comp_liq_prof) {
    // delta perturbation function for the liquid phase
    U.is(1,vl_indx,vl_indxe);
    v_l.set(U);
    U.rs();
    tau=double(tau_supg_c.get(vl_indx,vl_indx));
    tmp9.prod(v_l,grad_N,-1,-1,1).scale(tau);
    tmp6.prod(tmp9,Id,1,2,3);
    P_supg.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe);
    P_supg.add(tmp6);
    P_supg.rs();
    }

    if (comp_gas_prof) {
    // delta perturbation function for the gas phase
    U.is(1,vg_indx,vg_indxe);
    v_g.set(U);
    U.rs();
    tau=double(tau_supg_c.get(vg_indx,vg_indx));
    tmp9.prod(v_g,grad_N,-1,-1,1).scale(tau);
    tmp6.prod(tmp9,Id,1,2,3);
    P_supg.is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe);
    P_supg.add(tmp6);
    P_supg.rs();
    }

    if (comp_liq_prof) {
    // epsilon perturbation function for the liquid phase
    tau=double(tau_supg_c.get(vl_indx,vl_indx));
    rho_m=arho_l+arho_g;
    rho_m=rho_l;
    tau=tau/rho_m;
    tmp7.set(grad_N).scale(tau);
    tmp7.t();
    P_supg.ir(2,1).is(3,vl_indx,vl_indxe).add(tmp7);
    P_supg.rs();
    }

    if (comp_gas_prof) {

    if (0) {
    // epsilon perturbation function for the gas phase
    tau=double(tau_supg_c.get(2,2));
    rho_m=rho_l;
    tau=tau/rho_m;
    tmp7.rs();
    // Debug & fixme
    tmp7.set(grad_N).scale(tau);
    // end debug
    tmp7.t();
    P_supg.ir(2,2).is(3,vg_indx,vg_indxe).add(tmp7);
    P_supg.rs();
    } else {
    // new stabilization for gas continuity equation
    U.is(1,vg_indx,vg_indxe);
    v_g.set(U);
    U.rs();
    //    tau=double(tau_supg_c.get(vg_indx,vg_indx));
    tau=double(tau_supg_c.get(2,2));
    tmp9.prod(v_g,grad_N,-1,-1,1).scale(tau);

    // DEBUG :> fixme
    // tmp9.set(0.0);

    P_supg.ir(2,2).ir(3,2).add(tmp9).rs();
    //    P_supg.ir(2,2).ir(3,2).set(tmp9).rs();

    }
    }

    // added stabilization for continuity equation
    if (0) {
    // stabilization enhancement of continuity equation
    tau=double(tau_supg_c.get(1,1));
    rho_m=arho_l+arho_g;
    rho_m=rho_l;
    tau=tau*rho_m;
    tmp7.rs();
    tmp7.set(grad_N).scale(tau);
    tmp7.t();
    P_supg.is(2,vl_indx,vl_indxe).ir(3,1).add(tmp7);
    P_supg.rs();
    P_supg.is(2,vg_indx,vg_indxe).ir(3,2).add(tmp7);
    P_supg.rs();
    }
    tmp7.rs();
    // end debug
    }
}
#endif

//__INSERT_LICENSE__
//$Id: bubbly.cpp,v 1.16.8.1 2004/05/21 21:10:56 mstorti Exp $

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly.h"

// incluido
#include "./fm2funm.h"

class MyFun : public FastMat2_funm {
public:
  // No hardening at all
  double f(double l) { return sqrt(l); }
//  double f(double l) { return l; }
} my_fun;

// fin incluido

extern const char * jobinfo_fields;

extern GlobParam *GLOB_PARAM;

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

  //o Add LES for this particular elemset.
  EGETOPTDEF_ND(elemset,int,LES,0);
  //o van Driest constant for the damping law.
  EGETOPTDEF_ND(elemset,double,A_van_Driest,0);
  assert(A_van_Driest==0.);
  //o Smagorinsky constant.
  EGETOPTDEF_ND(elemset,double,C_smag,0.18);

  //o turbulence induced by bubble according to Sato model (constant)
  EGETOPTDEF_ND(elemset,double,Sato_model_coef,0.);

  //o stabilization factor (tau_fac)
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);
  //o shocap
  EGETOPTDEF_ND(elemset,double,shocap,1.);

  //o Adjust the stability parameters, taking into account
  // the time step. If the \verb+steady+ option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  // \verb+temporal_stability_factor+ is set to 0.
  EGETOPTDEF_ND(elemset,double,temporal_stability_factor,1.);
  // bubble diameter
  EGETOPTDEF_ND(elemset,double,d_bubble,0);
  // interphase drag force term key
  EGETOPTDEF_ND(elemset,int,comp_interphase_terms,0);
  // virtual mass key
  EGETOPTDEF_ND(elemset,int,comp_virtual_mass,0);
  // virtual mass coefficient
  EGETOPTDEF_ND(elemset,double,C_vm,0.5);

  // lift coefficient (null by default)
  EGETOPTDEF_ND(elemset,double,C_lift,0.0);

  // flag masking for matrices (Cp,Ajac,Cjac,Djac)
  // 1111 todas prendidas
  EGETOPTDEF_ND(elemset,int,flag_debug,0);
  EGETOPTDEF_ND(elemset,int,mask_matrix,1111);

  // key to add stabilization for gas continuity equation as for incompressible
  // as suggested by a Finland researcher
  EGETOPTDEF_ND(elemset,int,use_pspg_for_gas,0);
  EGETOPTDEF_ND(elemset,double,factor_pspg_for_gas,1.0);

  // key for using PMM method
  EGETOPTDEF_ND(elemset,int,use_pmm_method,0);

  // sonic velocity for remove singularity in Cp matrix
  EGETOPTDEF_ND(elemset,double,cc,1.e+6);

  // sonic velocity for remove singularity in Cp matrix
  EGETOPTDEF_ND(elemset,double,factor_sonic_speed,0.0);

  // to activate or not the grad alpha_g in source term
  EGETOPTDEF_ND(elemset,double,flag_grad_alpha_source,1.0);

  // to activate or not the interfacial pressure term
  EGETOPTDEF_ND(elemset,double,flag_interfacial_pressure,0.0);

  // to scale the temporal term in the liquid continuity equation
  EGETOPTDEF_ND(elemset,double,factor_liq_mass_eq_mod,1.0);

  // key for using different drag model
  EGETOPTDEF_ND(elemset,int,drag_model,0);

  // key for using different drag model
  EGETOPTDEF_ND(elemset,double,drag_value,0.);

  // key for using different drag model
  EGETOPTDEF_ND(elemset,int,upwind_gas_matrix,0);

  // to activate or not the turbulent dispersion term
  EGETOPTDEF_ND(elemset,double,coef_turbulent_dispersion,0.0);

  // Turbulent Schmidt number
  EGETOPTDEF_ND(elemset,double,Sc_t,1.0);

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

  grad_v_l_old.resize(2,ndim,ndim);
  grad_v_g_old.resize(2,ndim,ndim);

  grad_k.resize(1,ndim);
  grad_e.resize(1,ndim);
  grad_alpha_g.resize(1,ndim);
  grad_alpha_g_old.resize(1,ndim);
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

  if (comp_virtual_mass>0) {
//	  rho_l_eff= rho_l*(1.0-C_vm);
	  rho_l_eff= rho_l;
	  rho_g_eff= rho_g + C_vm*rho_l;
	  Dmat_vl.resize(1,ndim);
	  Dmat_vg.resize(1,ndim);
	  tmp_vm.resize(2,ndim,ndim);
	  tmp2_vm.resize(1,ndim);
  } else {
	  rho_l_eff = rho_l;          ;
	  rho_g_eff = rho_g;
  }


  if (use_pmm_method) {
    //    assert(!strcmp(jobinfo_fields,"gasliq"));
    Tau_beta.resize(2,ndof,ndof);
    my_fun.init(Tau_beta);
    Cpi.resize(2,ndof,ndof);
    // no se puede usar presion interfacial y pmm al mismo tiempo
    assert(flag_interfacial_pressure==0);
  }

  if (upwind_gas_matrix) {
    assert(ndim==2);
    sign_Lambda_gas.resize(2,ndim+1,ndim+1);
    VV_gas.resize(2,ndim+1,ndim+1);
    VVi_gas.resize(2,ndim+1,ndim+1);
    temp_tau_A.resize(2,ndim+1,ndim+1);
    temp2_tau_A.resize(2,ndim+1,ndim+1);
//    tau_A.resize(3,ndim,ndof,ndof);
    tau_A.resize(2,ndof,ndof);
    }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::element_hook(ElementIterator &element_a) {
	element=element_a;
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

#if 1
  // limito las fracciones de vacio entre 0 y 1
  double tol=1e-6;
  double alpha_g_ctf = (alpha_g < tol ? tol : alpha_g);
  alpha_g_ctf = (alpha_g_ctf > (1.0-tol) ? (1.0-tol) : alpha_g_ctf);
  double alpha_l_ctf = 1.0-alpha_g_ctf;
  alpha_g = alpha_g_ctf;
  alpha_l = alpha_l_ctf;
#endif

  arho_l = alpha_l*rho_l;
  arho_g = alpha_g*rho_g;

  arho_l_eff = alpha_l*rho_l_eff;
  arho_g_eff = alpha_g*rho_g_eff;

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

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");


  H.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  H.setel(arho_l,1);
  H.is(1,vl_indx,vl_indxe).set(v_l).scale(arho_l_eff);
  H.rs();
  }
  if (comp_kep_prof) {
  // k-eps turbulence
  H.setel(arho_l*k,k_indx);
  H.setel(arho_l*eps,e_indx);
  }
  if (comp_gas_prof) {
  // Gas phase
  if (use_pmm_method) {
	  compute_beta_pmm();
	  H.setel(arho_g*beta_pmm,2);
  } else {
      H.setel(arho_g,2);
  }
  H.is(1,vg_indx,vg_indxe).set(v_g).scale(arho_g_eff);
  H.rs();
  }

//  if (flag_debug) printf(" beta_pmm (comp enthalpy)= %15.5e \n",beta_pmm);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
	      double weight) {
  W_N.prod(W,N,1,2).scale(weight);
  W_Cp_N.prod(W_N,Cp,1,3,2,4);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg) {
  P_Cp.prod(P_supg,Cp,1,-1,-1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::get_C(FastMat2 &C) {
  C.set(Cjac);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::get_Cp(FastMat2 &Cpp) {
  Cpp.set(Cp);
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_tau_beta() {

  //    static FastMat2 svec,tmp9;

#define SHV(pp) pp.print(#pp ": ")

  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

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

  /*
    double Peclet = velmod * h_supg / (2. * visco_supg);
    double fz = (Peclet < 3. ? Peclet/3. : 1.);
    delta_supg = 0.5*h_supg*velmod*fz;
  */

//  compute_beta_pmm();

  SHV(U);

  tmp_beta_40.set(Cp);
  tmp_beta_40.setel(alpha_l/square(cc),1,1);
  tmp_beta_40.setel(beta_pmm*alpha_g/square(cc),2,1);
  if (comp_kep_prof==0) {
  tmp_beta_40.setel(1.,k_indx,k_indx);
  tmp_beta_40.setel(1.,e_indx,e_indx);
  }
  SHV(tmp_beta_40);
  Cpi.inv(tmp_beta_40);
  SHV(Cpi);

  Tau_beta.set(0.);
  for (int k=1; k<=ndim; k++) {
    Ajac.ir(1,k);
    SHV(Ajac);
    tmp_beta_41.prod(Cpi,Ajac,1,-1,-1,2);
    tmp_beta_42.prod(tmp_beta_41,tmp_beta_41,1,-1,-1,2);
    Tau_beta.add(tmp_beta_42);
  }
  Ajac.rs();

  if (comp_kep_prof==0) {
    Tau_beta.setel(1.,k_indx,k_indx);
    Tau_beta.setel(1.,e_indx,e_indx);
  }

  SHV(Tau_beta);

  my_fun.apply_ns(Tau_beta,fTau_beta);

  SHV(fTau_beta);

  tmp_beta_4.inv(fTau_beta).scale(0.5*h_supg);


  Tau_beta.prod(tmp_beta_4,Cpi,1,-1,-1,2);

  SHV(Tau_beta);

  /*
  PetscPrintf(PETSC_COMM_WORLD," ENDING \n");
  PetscFinalize();
  exit(0);
  */

#undef SHV

}

/*
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_beta_pmm() {

  double tol=1.0e-10;

  v_rel.set(v_l).rest(v_g);
  tmp_beta_1.prod(v_g,v_g,-1,-1);
  double v1 = double(tmp_beta_1)*rho_g*alpha_l;
  tmp_beta_2.set(v_g).scale(rho_g*alpha_l);
  tmp_beta_3.set(v_l).scale(rho_l*alpha_g);
  tmp_beta_2.rest(tmp_beta_3);

  double v2 = -(double(tmp_beta_31.prod(tmp_beta_2,v_g,-1,-1)) +
		double(tmp_beta_32.prod(v_rel,v_l,-1,-1).scale(rho_l*alpha_g)));
  double v3 = -double(tmp_beta_1)*rho_l*alpha_g;

  double v_rel_mod = sqrt(v_rel.sum_square_all());
  FastMat2::branch();
  if ((v_rel_mod<tol) || (v1+v2+v3)>0) {
    FastMat2::choose(0);
    beta_pmm = 1.0;
  } else {
    double disc = sqrt(square(v2)-4.0*v1*v3);
    double beta_1 = (-v2+disc)/(2.*v1);
    double beta_2 = (-v2-disc)/(2.*v1);
    double beta_max = beta_1; if(beta_2>beta_1) beta_max = beta_2;
    double beta_min = beta_1; if(beta_2<beta_1) beta_min = beta_2;
    double factor = 0.01;
    //    FastMat2::branch();
    if (fabs(beta_min-1.0) < fabs(beta_max-1.0)) {
      //      FastMat2::choose(0);
      beta_pmm = beta_min - factor*fabs(beta_min);
    } else {
      beta_pmm = beta_max + factor*fabs(beta_max);
    }
    //    FastMat2::leave();
  }
  FastMat2::leave();

  // DEBUG
  // printf(" beta_pmm (comp beta)= %15.5e \n",beta_pmm);
  // beta_pmm = 2.0;
  // END DEBUG
}
*/

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_beta_pmm() {

//  double v1 = double(tmp_beta_1.prod(v_g,v_g,-1,-1));
//  v1 = sqrt(v1);

//  FastMat2 ugn,uln;

  double tol=1.0e-10;

  v_rel.set(v_l).rest(v_g);
  double v_rel_mod = sqrt(v_rel.sum_square_all());

  double v1 = sqrt(v_g.sum_square_all());
  tmp_beta_2.set(v_g).scale(1.0/v1);
  ugn.prod(v_g,tmp_beta_2,-1,-1);
  uln.prod(v_l,tmp_beta_2,-1,-1);

  double deno = rho_g*alpha_l*square(ugn);

  double beta_1,beta_2;

  FastMat2::branch();
  if (deno>tol) {
    FastMat2::choose(0);
  double xi = 1;
  double nume = rho_g*alpha_l*square(ugn)+2*rho_l*alpha_g*square(uln)-2*rho_l*alpha_g*uln*ugn+
        xi*2*sqrt(rho_l*alpha_g*(square(ugn)*rho_g*alpha_l+square(uln)*rho_l*alpha_g)*square(uln-ugn));
  beta_1 = nume/deno;

  xi = -1;
  nume = rho_g*alpha_l*square(ugn)+2*rho_l*alpha_g*square(uln)-2*rho_l*alpha_g*uln*ugn+
        xi*2*sqrt(rho_l*alpha_g*(square(ugn)*rho_g*alpha_l+square(uln)*rho_l*alpha_g)*square(uln-ugn));
  beta_2 = nume/deno;
} else {
	beta_1 = 1;
	beta_2 = 1;
}
  FastMat2::leave();

  FastMat2::branch();
  if (v_rel_mod<tol) {
    FastMat2::choose(0);
    beta_pmm = 1.0;
  } else {
    double beta_max = beta_1; if(beta_2>beta_1) beta_max = beta_2;
    double beta_min = beta_1; if(beta_2<beta_1) beta_min = beta_2;
    double factor = 0.01;
    if (fabs(beta_min-1.0) < fabs(beta_max-1.0)) {
      beta_pmm = beta_min - factor*fabs(beta_min);
    } else {
      beta_pmm = beta_max + factor*fabs(beta_max);
    }
  }
  FastMat2::leave();

  beta_pmm = (fabs(beta_pmm) < 1.0e-6 ? 1.0e-6 : beta_pmm);

  // DEBUG
  // printf(" beta_pmm (comp beta)= %15.5e \n",beta_pmm);
  // beta_pmm = 2.0;
  // END DEBUG
}

/*
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
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

    // DEBUG
    if (0) {
      h_supg = 0.4/100.;
      h_pspg = h_supg;
    }


    //SHV(h_pspg);
    //SHV(h_supg);
    //SHV(velmod);

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

*/


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_tau(int ijob) {

  //    static FastMat2 svec,tmp9;

    const FastMat2 &grad_N = *advdf_e->grad_N();

    double tol=1.0e-12;
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


    double velmax = velmod; // + factor_sonic_speed*sqrt(3.0)*sqrt(cc_supg);
    double v1 = sqrt(v_g.sum_square_all());

//    double vel_bubbly_wave ;
    FastMat2::branch();
    if (ijob==1 && v1>tol) {
        FastMat2::choose(0);
        tmp_beta_2.set(v_g).scale(1.0/v1);
        ugn.prod(v_g,tmp_beta_2,-1,-1);
        uln.prod(v_l,tmp_beta_2,-1,-1);

		double nume = rho_g*alpha_l*square(ugn)+rho_l*alpha_g*square(uln);
		double deno = 0.5*(rho_g*alpha_l*ugn*(1+beta_pmm)+2*rho_l*alpha_g*uln);

		if (fabs(deno)>tol) {
		double vel_bubbly_wave = fabs(nume/deno);
		velmax = (vel_bubbly_wave > velmax ? vel_bubbly_wave : velmax);
	    }
    if (flag_debug) U.print("U:\n");
    if (flag_debug) printf(" compute_tau = %15.7e %15.7e %15.7e \n",beta_pmm,velmod,velmax);
//    if (flag_debug) printf(" compute_tau = %15.7e %15.7e %15.7e \n",alpha_g,double(ugn),double(uln));

	}
    FastMat2::leave();

    FastMat2::branch();
	if (ijob==2 && v1>tol) {
        FastMat2::choose(0);
        tmp_beta_2.set(v_g).scale(1.0/v1);
        ugn.prod(v_g,tmp_beta_2,-1,-1);
        uln.prod(v_l,tmp_beta_2,-1,-1);
        double vr2 =square(uln-ugn);
        double disc = rho_g*vr2*(alpha_g*(rho_l*alpha_g-rho_g)+rho_g);
		disc = (disc > 0.0 ? disc : 0.0);
        double nume = rho_g*alpha_l*square(ugn)+rho_l*alpha_g*square(uln)-rho_g*vr2;
        double deno (rho_g*alpha_l*ugn+rho_l*alpha_g*uln+sqrt(disc));
        double L1,L2;
        if (fabs(deno)>tol) {
			L1 = fabs(nume/deno);
		} else {
			L1 = 0;
		}
        deno = (rho_g*alpha_l*ugn+rho_l*alpha_g*uln-sqrt(disc));
        if (fabs(deno)>tol) {
			L2 = fabs(nume/deno);
		} else {
			L2 = 0;
		}
        if (L2>L1)L1=L2;
		if (L1>0) {
		velmax = (L1 > velmax ? L1 : velmax);
	    }
    if (flag_debug) U.print("U:\n");
    if (flag_debug) printf(" compute_tau = %15.7e %15.7e %15.7e \n",L1,velmod,velmax);

    }
    FastMat2::leave();

    // DEBUG

    // tau_supg_a =  h_supg/2./velmax;
    // tau_supg_a =  h_supg/2./velmod;
    double velmod2=velmod;
    if (ijob==1) velmod2=velmod/fabs(beta_pmm);
    velmod2 = (velmod2 > velmod ? velmod2 : velmod);
    tau_supg_a =  square(2.*velmod2/h_supg)
                  + 9.*square(4.*visco_supg/square(h_supg));
    tau_supg_a = 1./sqrt(tau_supg_a);

    double pspg_factor=1.;
    tau_pspg = pspg_factor*tau_supg_a;

    double Peclet = velmax * h_supg / (2 *visco_supg);
    double fz = (Peclet < 3. ? Peclet/3. : 1.);
    delta_supg = 0.5*h_supg*velmax*fz;

	}


/*
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_tau_gas() {

  //    static FastMat2 svec,tmp9;

    double Cpi = flag_interfacial_pressure;
    double velmax = velmod; // + factor_sonic_speed*sqrt(3.0)*sqrt(cc_supg);
    double v1 = sqrt(vel_supg.sum_square_all());
    tmp_beta_2.set(vel_supg).scale(1.0/v1);
    ugn.prod(v_g,tmp_beta_2,-1,-1);
    uln.prod(v_l,tmp_beta_2,-1,-1);
    double vr2 = square(uln-ugn);
    double tau;

// EN CONSTRUCCION
    double vel_wave=sqrt(Cpi*vr2);


//    printf(" vel wave  = %15.5e \n",vel_wave);

    double vel;

    tau_A.set(0.0);

    FastMat2::branch();
    if (vel_wave>1.0e-16) {
      FastMat2::choose(0);
      for (int j=1; j<=ndim; j++) {
        vel=double(vel_supg.get(j));
        sign_Lambda_gas.set(0.0);
        double vaux = (fabs(vel) < 1.0e-16 ? 0 : vel/fabs(vel) );
        sign_Lambda_gas.setel(vaux,1,1);
        vaux = (fabs(vel+vel_wave) < 1.0e-16 ? 0 : (vel+vel_wave)/fabs(vel+vel_wave) );
        sign_Lambda_gas.setel(vaux,2,2);
        vaux = (fabs(vel-vel_wave) < 1.0e-16 ? 0 : (vel-vel_wave)/fabs(vel-vel_wave) );
        sign_Lambda_gas.setel(vaux,3,3);

	VV_gas.set(0.0);
	VV_gas.setel(1.,1,2);
	VV_gas.setel(1.,1,3);
	VV_gas.setel(1.,3,1);
	VV_gas.setel(vel_wave/alpha_g,2,2);
	VV_gas.setel(-vel_wave/alpha_g,2,3);

	VVi_gas.set(0.0);
	VVi_gas.setel(1.,1,3);
	VVi_gas.setel(0.5,2,1);
	VVi_gas.setel(0.5,3,1);
	VVi_gas.setel(0.5*alpha_g/vel_wave,2,2);
	VVi_gas.setel(-0.5*alpha_g/vel_wave,3,2);

	// masking tau_supg for all gas variables

	tau_A.ir(1,j).is(2,2).is(2,vg_indx,vg_indxe).is(3,2).is(3,vg_indx,vg_indxe);
	temp_tau_A.prod(VV_gas,sign_Lambda_gas,1,-1,-1,2);
	temp2_tau_A.prod(temp_tau_A,VVi_gas,1,-1,-1,2).scale(0.5*h_supg);
	tau_A.add(temp2_tau_A).rs();
      }

    } else {

      tau=double(tau_supg_c.get(vg_indx,vg_indx));

      tau_A.is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe);
      tau_A.prod(v_g,Id,1,2,3).scale(tau).rs();

      tau_A.ir(2,2).ir(3,2);
      tau_A.add(v_g).scale(tau).rs();

    }

    FastMat2::leave();
    tau_A.scale(tau_fac);

}

*/

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::compute_tau_gas() {


    double tol=1.0e-10;

    double Cpi = flag_interfacial_pressure;

    double rho_vm = 1.0;
    if (comp_virtual_mass>0) rho_vm = rho_g/(C_vm*rho_l+rho_g);

    double velmax = velmod; // + factor_sonic_speed*sqrt(3.0)*sqrt(cc_supg);
    double v1 = sqrt(vel_supg.sum_square_all());
    tmp_beta_2.set(vel_supg).scale(1.0/v1);
    ugn.prod(v_g,tmp_beta_2,-1,-1);
    uln.prod(v_l,tmp_beta_2,-1,-1);
    double vr2 = square(uln-ugn);
    double tau;

// EN CONSTRUCCION
    double vel_wave=sqrt(Cpi*vr2*rho_vm);


//    printf(" vel wave  = %15.5e \n",vel_wave);

    double vaux ,vel ,nx ,ny;

    nx = double(tmp_beta_2.get(1));
    ny = double(tmp_beta_2.get(2));

// en realidad tau_A ===> para esta version es solo la matriz tau
//             sign_Lambda_gas ===> son las inversas de los valores absolutos de los autovalores

    tau_A.set(0.0);

    FastMat2::branch();
    if (vel_wave>tol && v1> tol) {
      FastMat2::choose(0);
//      for (int j=1; j<=ndim; j++) {
//        vel=double(vel_supg.get(j));
        sign_Lambda_gas.set(0.0);
        vaux = (fabs(v1+vel_wave) < tol ? tol : fabs(v1+vel_wave) );
        sign_Lambda_gas.setel(1./vaux,1,1);
        vaux = (fabs(v1-vel_wave) < tol ? tol : fabs(v1-vel_wave) );
        sign_Lambda_gas.setel(1./vaux,2,2);
        double vaux = (fabs(v1) < tol ? tol : fabs(v1) );
        sign_Lambda_gas.setel(1./vaux,3,3);

	VV_gas.set(0.0);
	VV_gas.setel(1.,1,1);
	VV_gas.setel(1.,1,2);

	VV_gas.setel(vel_wave/alpha_g*nx,2,1);
	VV_gas.setel(-vel_wave/alpha_g*nx,2,2);
	VV_gas.setel(-ny,2,3);

	VV_gas.setel(vel_wave/alpha_g*ny,3,1);
	VV_gas.setel(-vel_wave/alpha_g*ny,3,2);
	VV_gas.setel(nx,3,3);

	VVi_gas.set(0.0);
	VVi_gas.setel(0.5,1,1);
	VVi_gas.setel(0.5*alpha_g/vel_wave*nx,1,2);
	VVi_gas.setel(0.5*alpha_g/vel_wave*ny,1,3);

	VVi_gas.setel( 0.5,2,1);
	VVi_gas.setel(-0.5*alpha_g/vel_wave*nx,2,2);
	VVi_gas.setel(-0.5*alpha_g/vel_wave*ny,2,3);

	VVi_gas.setel(-ny,3,2);
	VVi_gas.setel( nx,3,3);

	// masking tau_supg for all gas variables

	tau_A.is(1,2).is(1,vg_indx,vg_indxe).is(2,2).is(2,vg_indx,vg_indxe);
	temp_tau_A.prod(VV_gas,sign_Lambda_gas,1,-1,-1,2);
	temp2_tau_A.prod(temp_tau_A,VVi_gas,1,-1,-1,2).scale(0.5*h_supg);
	tau_A.add(temp2_tau_A).rs();

    }

    FastMat2::leave();
    tau_A.scale(tau_fac);

/*
    U.print("U:\n");
    VV_gas.print("V:\n");
    VVi_gas.print("Vi:\n");
    Ajac.print("Ajac:\n");
    sign_Lambda_gas.print("inv(abs(La)):\n");

    PetscPrintf(PETSC_COMM_WORLD," ENDING \n");
    PetscFinalize();
    exit(0);
*/

}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
// compute the profile for each specific subproblem according to jobinfo value
void bubbly_ff::set_profile(FastMat2 &seed) {

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  //  Matrix seed;
  //  seed= Matrix(ndof,ndof);
  seed.set(0.);

  if (comp_liq_prof && comp_gas_prof) {
    for (int j=1; j<=2*ndim+2; j++) {
      for (int k=1; k<=2*ndim+2; k++) {
        seed.setel(1.,j,k);
      }
    }
  }
  else if (comp_liq_prof) {
    int ip[] = {1,3,4,5};
    for (int j=0; j<=ndim; j++) {
      for (int k=0; k<=ndim; k++) {
	//        seed(ip[j],ip[k])=1;
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  }

  else if (comp_gas_prof) {
    int ip[] = {2,2+ndim+1,2+ndim+2,2+ndim+3};
    for (int j=0; j<=ndim; j++) {
      for (int k=0; k<=ndim; k++) {
	//        seed(ip[j],ip[k])=1;
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  }

  else if (comp_kep_prof) {
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

  double coe,pot,f_drag_alpha , df_drag_alpha ;

  options &= ~SCALAR_TAU;	// tell the advective element routine
				// that we are returning a MATRIX tau

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  double Dt = GLOB_PARAM->Dt;
  double rec_Dt = 1./Dt;
  if (GLOB_PARAM->steady) rec_Dt = 0.;

  int mask_Cp=1, mask_Ajac=1, mask_Cjac=1, mask_Djac=1;

  if(flag_debug){
  // mask for each matrix in order to activate or not for debugging purposes
  mask_Cp = mask_matrix/1000;
  mask_Ajac = (mask_matrix-1000*mask_Cp)/100;
  mask_Cjac = (mask_matrix-1000*mask_Cp-100*mask_Ajac)/10;
  mask_Djac = (mask_matrix-1000*mask_Cp-100*mask_Ajac-10*mask_Cjac);
}

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Enthalpy Jacobian Cp

  Cp.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Cp.setel(-rho_l*factor_liq_mass_eq_mod,1,2);
//  Cp.is(1,vl_indx,vl_indxe).ir(2,2).set(v_l).scale(-rho_l_eff).rs();
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(arho_l_eff).rs();
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
  Cp.setel(rho_g,2,2);
//  Cp.is(1,vg_indx,vg_indxe).ir(2,2).set(v_g).scale(rho_g_eff).rs();
  Cp.rs().is(1,vg_indx,vg_indxe).is(2,vg_indx,vg_indxe).eye(arho_g_eff).rs();

  // Cp.print("Cp:\n");
  if (use_pmm_method) {
	  compute_beta_pmm();
      Cp.ir(1,2).scale(beta_pmm).rs();
//      if (flag_debug) printf(" beta_pmm (comp Cp)= %15.5e \n",beta_pmm);
  }
  }

  // Cp.print("Cp:\n");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Advective fluxes
  flux.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  flux.ir(1,1).set(v_l).scale(arho_l).rs();
  // Amoml.prod(v_l,v_l,1,2).scale(rho_l_eff).axpy(Id,p);
  Amoml.prod(v_l,v_l,1,2).scale(rho_l_eff);
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
  //  Amomg.prod(v_g,v_g,1,2).scale(rho_g);
  // Amomg.prod(v_g,v_g,1,2).scale(rho_g_eff).axpy(Id,p);
  Amomg.prod(v_g,v_g,1,2).scale(rho_g_eff);
  flux.is(1,vg_indx,vg_indxe).axpy(Amomg,alpha_g).rs();
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Adjective Jacobians
  Ajac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Ajac.rs().ir(2,1).ir(3,2).set(v_l).scale(-rho_l);
  Ajac.rs().ir(2,1).is(3,vl_indx,vl_indxe).axpy(Id,arho_l);

  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).axpy(Id,alpha_l);

  if (0) {
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,2).axpy(Amoml,-1);
  Y.rs().prod(v_l,Id,2,1,3).scale(arho_l_eff);
  Ajac.rs().is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe)
    .add(Y);
  Y.prod(v_l,Id,1,2,3).scale(arho_l_eff).rs();
  Ajac.add(Y);
  Y.rs();
  Ajac.rs();
} else {
  // escrito para usarlo con la forma no conservativa
  // Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,2).axpy(Id,-p);
  Y.prod(v_l,Id,1,2,3).scale(arho_l_eff).rs();
  Ajac.rs().is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).add(Y);
  Ajac.rs();
}
  }

  if (comp_gas_prof) {
  // Gas phase
  Ajac.rs().ir(2,2).ir(3,2).set(v_g).scale(rho_g);
  Ajac.rs().ir(2,2).is(3,vg_indx,vg_indxe).axpy(Id,arho_g);

//  Ajac.rs().is(2,vg_indx,vg_indxe).ir(3,1).axpy(Id,alpha_g);

  if (0) {
  Ajac.rs().is(2,vg_indx,vg_indxe).ir(3,2).set(Amomg);
  Y.prod(v_g,Id,2,1,3).scale(arho_g_eff);
  Ajac.rs().is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe)
    .add(Y);
  Y.prod(v_g,Id,1,2,3).scale(arho_g_eff).rs();
  Ajac.add(Y);
  Y.rs();
  Ajac.rs();
  } else {
  Y.prod(v_g,Id,1,2,3).scale(arho_g_eff).rs();
  Ajac.rs().is(2,vg_indx,vg_indxe).is(3,vg_indx,vg_indxe).add(Y);
  Ajac.rs();
  }
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

  // Ajac.print("Ajac:\n");


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

if (0) {
  // Turbulent viscosity
  // fixme:= la tasa de deformacion se toma la del liquido??
  strain_rate_scalar = strain_rate_l.sum_square_all();
  // fixme:= esto no compila
  visco_t = C_mu* square(k)/eps;
  P_k = 2*rho_l*visco_t*strain_rate_scalar;
}
	// Smagorinsky turbulence model
//	double nu_eff;
	if ((LES) && (options & COMP_UPWIND)) {

    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);
#define pi M_PI
    double Volume = advdf_e->volume();
    double Delta;
    if (ndim==2) Delta = sqrt(Volume);
    if (ndim==3) Delta = cbrt(Volume);

	  double tr = (double) tmp15.prod(strain_rate_l,strain_rate_l,-1,-2,-1,-2);
	  double van_D;
/*
	  if (A_van_Driest>0.) {
	    dist_to_wall.prod(SHAPE,xloc,-1,-1,1).rest(wall_coords);
	    double ywall = sqrt(dist_to_wall.sum_square_all());
	    double y_plus = ywall*shear_vel/VISC;
	    van_D = 1.-exp(-y_plus/A_van_Driest);
	    if (k % 250==0) printf("van_D: %f\n",van_D);
	  } else van_D = 1.;
*/
      van_D = 1.;

	  visco_t = SQ(C_smag*Delta*van_D)*sqrt(2*tr);
	} else {
	  visco_t = 0.;
	}

        double visco_sato = 0.;
	if (Sato_model_coef>0) {
          visco_sato = Sato_model_coef*rho_l*alpha_g*d_bubble*vslip;
	}

  // limito la viscosidad por debajo para evitar valores negativos
  double tol=0.0;
  double alpha_g_ctf = (alpha_g < tol ? tol : alpha_g);
  alpha_g_ctf = (alpha_g_ctf > (1.0-tol) ? (1.0-tol) : alpha_g_ctf);
  double alpha_l_ctf = 1.0-alpha_g_ctf;

  visco_l_eff = visco_l + visco_sato + rho_l*visco_t;
  visco_g_eff = visco_g + rho_g*visco_t + rho_g/rho_l*visco_sato;

  // Strain rate for the gas
  grad_U.is(2,vg_indx,vg_indxe);
  grad_v_g.set(grad_U);
  grad_U.rs();
  strain_rate_g.set(grad_v_g);
  grad_v_g.t();
  strain_rate_g.add(grad_v_g).scale(0.5);
  grad_v_g.rs();

  // gas velocity divergence
  grad_v_g.d(1,2);
  tmp10.sum(grad_v_g,-1);
  grad_v_g.rs();
  double div_vg = double(tmp10);

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

  grad_U.ir(2,2);
  grad_alpha_g.set(grad_U);
  grad_U.rs();
  grad_alpha_l.set(grad_alpha_g).scale(-1.0);

  grad_U.ir(2,1);
  grad_p.set(grad_U);
  grad_U.rs();

  v_g_l.set(v_g).rest(v_l);
  v_slip = v_g_l.sum_square_all();
  v_slip = sqrt(v_slip);

  // Reactive terms
  //   [1] bouyancy terms
  //   [2] Phase interaction forces
  //   [3] pressure gradient terms for gas phase
  Cjac.set(0.);
  G_source.set(0.);

  if ((options & COMP_SOURCE) || (options & COMP_SOURCE_LUMPED)){

  if (comp_interphase_terms==1) {
  assert(d_bubble>0);
  C1_drag = 3./4.*visco_l/d_bubble/d_bubble;
  Rey_bubble =rho_l*d_bubble*v_slip/visco_l;

    FastMat2::branch();
    if (Rey_bubble>0){
      FastMat2::choose(0);

	  // user defined C_drag function (a pata por ahora)
      // Lai & Salcudean
      if (drag_model==0) {
      C_drag_ff= 24.0/Rey_bubble;
      dCDdRe_ff= -24.0/Rey_bubble/Rey_bubble;
      } else if (drag_model==1) {
      double tmp0_drag = (1.0+0.15*pow(Rey_bubble,0.687));
	  C_drag_ff= 24.0/Rey_bubble*tmp0_drag;
	  dCDdRe_ff= -24.0/Rey_bubble/Rey_bubble*tmp0_drag +
                      24.0/Rey_bubble*(0.15*0.687*pow(Rey_bubble,0.687-1.0));
  } else if (drag_model==2) {
//	  C_drag_ff= 0.44;
	  C_drag_ff= drag_value;
	  dCDdRe_ff= 0.0;
  }

	  dRedU = rho_l*d_bubble/visco_l;

	  tmp1_drag = C1_drag*Rey_bubble*C_drag_ff;

      coe = 1.0;
      pot = 1.0;
      f_drag_alpha = coe*pow(alpha_l,pot)*alpha_g;
      df_drag_alpha = coe*pow(alpha_l,pot-1)*(1.-(1.+pot)*alpha_g);

      tmp2_drag = tmp1_drag*df_drag_alpha;

          tmp4_drag = C_drag_ff+Rey_bubble*dCDdRe_ff;

/*
          Phi_1.prod(v_g_l,v_g_l,1,2).scale(-tmp4_drag*dRedU/v_slip);
          Phi_2.set(Id).scale(Rey_bubble*C_drag_ff);
          Phi_1.rest(Phi_2);
*/
          Phi_1.prod(v_g_l,v_g_l,1,2).scale(tmp4_drag*dRedU/v_slip);
          Phi_2.set(Id).scale(Rey_bubble*C_drag_ff);
          Phi_1.add(Phi_2);

	  //Phi_1.scale(C1_drag*alpha_l*alpha_g);
	  //Phi_1.scale(C1_drag*alpha_g);
	  Phi_1.scale(C1_drag*f_drag_alpha);

	tmp3_drag = tmp1_drag*f_drag_alpha;


          if (comp_liq_prof) {
          // Liquid phase
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).axpy(Phi_1,id_liq*id_liq).rs();
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vg_indx,vg_indxe).axpy(Phi_1,id_liq*id_gas).rs();
	  Cjac.is(1,vl_indx,vl_indxe).ir(2,2).set(v_g_l).scale(-tmp2_drag*id_liq).rs();
          }
          if (comp_gas_prof) {
          // Gas phase
	  Cjac.is(1,vg_indx,vg_indxe).is(2,vg_indx,vg_indxe).axpy(Phi_1,id_gas*id_gas).rs();
	  Cjac.is(1,vg_indx,vg_indxe).is(2,vl_indx,vl_indxe).axpy(Phi_1,id_gas*id_liq).rs();
	  Cjac.is(1,vg_indx,vg_indxe).ir(2,2).set(v_g_l).scale(-tmp2_drag*id_gas).rs();
	  }


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

    // LIFT
    tmp1_lift.set(v_g_l).scale(C_lift*rho_l*alpha_g);
    rotor_v_l.set(grad_v_l);
    grad_v_l.t();
    rotor_v_l.rest(grad_v_l);
    grad_v_l.rs();
    tmp2_lift.prod(tmp1_lift,rotor_v_l,-1,1,-1);

    if (flag_debug) tmp2_lift.print("Lift  :");

        if (comp_liq_prof) {
        // Liquid phase
	G_source.is(1,vl_indx,vl_indxe).add(tmp2_lift).rs();
	}
        if (comp_gas_prof) {
        // Gas phase
	G_source.is(1,vg_indx,vg_indxe).axpy(tmp2_lift,-1.0).rs();
	}

// turbulent dispersion
    double CTD = coef_turbulent_dispersion*tmp1_drag*visco_t/Sc_t;

    if (comp_liq_prof) {
//      G_source.is(1,vl_indx,vl_indxe).axpy(grad_alpha_g,-CTD).rs();
      Ajac.is(2,vl_indx,vl_indxe).ir(3,2).axpy(Id,CTD).rs();
      }

    if (comp_gas_prof) {
//      G_source.is(1,vg_indx,vg_indxe).axpy(grad_alpha_g,CTD).rs();
      Ajac.is(2,vg_indx,vg_indxe).ir(3,2).axpy(Id,-CTD).rs();
      }

  }
}


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

    double alpha_g_old = Uo.get(2);

    const FastMat2 &grad_N = *advdf_e->grad_N();

    //    double tau_supg_a,tau_pspg,delta_supg,visco_supg,velmod, h_supg;
    // Fase liquida
    vel_supg.set(v_l_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }

    visco_supg = visco_l_eff/rho_l;
    // cc_supg    = p/rho_l;
    // if (flag_debug) printf(" cc_supg (liq) = %15.5e \n",cc_supg);
    int ijob =0;
    compute_tau(ijob);

    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }
    tau_supg.eye(tau_supg_a).setel(delta_supg,1,1);

    // if (use_pmm_method) compute_tau_beta();
    // Fase gas
    vel_supg.set(v_g_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }
    visco_supg = visco_g_eff/rho_g;
    // cc_supg=p/rho_g;
    // if (flag_debug) printf(" cc_supg (gas) = %15.5e \n",cc_supg);
    ijob=0;
    if (use_pmm_method) {
		ijob=1;
    } else if (flag_interfacial_pressure>0) {
		ijob=2;
	}

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

    if (shocap>0) {

    double Dnum = double(tau_supg_c.get(2,2));;
    Dnum = Dnum*rho_g*shocap;

    fluxd.ir(1,2).set(grad_alpha_g).scale(Dnum);
    fluxd.rs();
    Djac.ir(2,2).ir(4,2).set(Id).scale(Dnum);
    Djac.rs();

  /*
    grad_U.is(2,vg_indx,vg_indxe);
    grad_v_g.set(grad_U);
    grad_U.rs();

    fluxd.is(1,vg_indx,vg_indxe).axpy(grad_v_g,alpha_g_old*Dnum);
    fluxd.rs();

    for (int k=vg_indx; k<=vg_indxe; k++) {
           Djac.ir(2,k).ir(4,k).axpy(Id,alpha_g_old*Dnum).rs();
    }
*/
    }
    }
  }

  if ((options & COMP_SOURCE) || (options & COMP_SOURCE_NOLUMPED)){
    // Bouyancy forces

    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);

    FastMat2 &Uo = (FastMat2 &) advdf_e->Uold();
    FastMat2 &grad_Uo = (FastMat2 &) advdf_e->grad_Uold();

    grad_Uo.ir(2,2);
    grad_alpha_g_old.set(grad_Uo);
    grad_Uo.rs();
    grad_alpha_l_old.set(grad_alpha_g_old).scale(-1.0);
    double pold = Uo.get(1);
    Uo.is(1,vl_indx,vl_indxe);
    v_l_old.set(Uo);
    Uo.rs();
    Uo.is(1,vg_indx,vg_indxe);
    v_g_old.set(Uo);
    Uo.rs();

    grad_Uo.is(2,vl_indx,vl_indxe);
    grad_v_l_old.set(grad_Uo);
    grad_Uo.rs();

    grad_Uo.is(2,vg_indx,vg_indxe);
    grad_v_g_old.set(grad_Uo);
    grad_Uo.rs();

    double alpha_g_old = Uo.get(2);
    double alpha_l_old = 1-alpha_g_old;

    if (comp_virtual_mass>0) {

//	  Dmat_vl.prod(v_l,grad_v_l,-1,1,-1).axpy(v_l,rec_Dt).axpy(v_l_old,-rec_Dt);
//	  Dmat_vg.prod(v_g,grad_v_g,-1,1,-1).axpy(v_g,rec_Dt).axpy(v_g_old,-rec_Dt);

	  Dmat_vl.prod(v_l_old,grad_v_l_old,-1,1,-1).axpy(v_l,rec_Dt).axpy(v_l_old,-rec_Dt);
	  Dmat_vg.prod(v_g_old,grad_v_g_old,-1,1,-1).axpy(v_g,rec_Dt).axpy(v_g_old,-rec_Dt);

    }

    double coef_interfacial_pressure = -rho_g*flag_interfacial_pressure*pow(v_slip,2.0);


//    if (flag_debug) printf(" coef_interfacial_pressure = %15.5e \n",coef_interfacial_pressure);
    if (comp_liq_prof) {
    // Liquid phase

//      G_source.is(1,vl_indx,vl_indxe).axpy(grad_alpha_g,-coef_interfacial_pressure).rs();
      Ajac.is(2,vl_indx,vl_indxe).ir(3,2).axpy(Id,coef_interfacial_pressure).rs();

     if (comp_virtual_mass>0) {
	  tmp2_vm.set(Dmat_vg).rest(Dmat_vl);

//	  tmp2_vm.scale(rho_l*C_vm*alpha_g);
	  tmp2_vm.scale(rho_l*C_vm*alpha_g_old);

	  G_source.is(1,vl_indx,vl_indxe).add(tmp2_vm).rs();
	  }
    }
    if (comp_gas_prof) {
    // Gas phase
      G_source.is(1,vg_indx,vg_indxe).axpy(grad_p,-alpha_g*flag_grad_alpha_source).rs();

//      G_source.is(1,vg_indx,vg_indxe).axpy(grad_alpha_g,coef_interfacial_pressure).rs();
      Ajac.is(2,vg_indx,vg_indxe).ir(3,2).axpy(Id,-coef_interfacial_pressure).rs();

     if (comp_virtual_mass>0) {
//	  tmp2_vm.set(Dmat_vl).scale(alpha_g*rho_l*C_vm);
  	  tmp2_vm.set(Dmat_vl).scale(alpha_g_old*rho_l*C_vm);
	  G_source.is(1,vg_indx,vg_indxe).add(tmp2_vm).rs();
	  }

    }


  if (comp_liq_prof) {
    //    Cjac.is(1,vl_indx,vl_indxe).ir(2,1).axpy(grad_alpha_l,-flag_grad_alpha_source).rs();


  if (comp_virtual_mass>0) {
    tmp_vm.set(grad_v_l).scale(0.0).axpy(Id,rec_Dt).scale(rho_l*C_vm*alpha_g_old);
    Cjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).add(tmp_vm).rs();
    tmp_vm.set(grad_v_g).scale(0.0).axpy(Id,rec_Dt).scale(-rho_l*C_vm*alpha_g_old);
    Cjac.is(1,vl_indx,vl_indxe).is(2,vg_indx,vg_indxe).add(tmp_vm).rs();
//    tmp2_vm.set(Dmat_vg).rest(Dmat_vl).scale(-rho_l*C_vm);
//    Cjac.is(1,vl_indx,vl_indxe).ir(2,2).add(tmp2_vm).rs();
	  }

  }


  if (comp_gas_prof) {
//    Cjac.is(1,vg_indx,vg_indxe).ir(2,2).axpy(grad_alpha_g,flag_grad_alpha_source).rs();
    Cjac.is(1,vg_indx,vg_indxe).ir(2,2).axpy(grad_p,flag_grad_alpha_source).rs();


  if (comp_virtual_mass>0) {
      tmp_vm.set(grad_v_l).scale(0.0).axpy(Id,rec_Dt).scale(-rho_l*C_vm*alpha_g_old);
      Cjac.is(1,vg_indx,vg_indxe).is(2,vl_indx,vl_indxe).add(tmp_vm).rs();
//    tmp2_vm.set(Dmat_vl).scale(-rho_l*C_vm);
//    Cjac.is(1,vg_indx,vg_indxe).ir(2,2).add(tmp2_vm).rs();
	  }

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


        if (comp_liq_prof) {
    Cjac.is(1,vl_indx,vl_indxe).ir(2,2).axpy(G_body,rho_l).rs();
    G_source.is(1,vl_indx,vl_indxe).axpy(G_body,arho_l).rs();
//    G_source.is(1,vl_indx,vl_indxe).axpy(G_body,rho_l*alpha_l_old).rs();
}
        if (comp_gas_prof) {
      Cjac.is(1,vg_indx,vg_indxe).ir(2,2).axpy(G_body,-rho_g).rs();
      G_source.is(1,vg_indx,vg_indxe).axpy(G_body,arho_g).rs();
//      G_source.is(1,vg_indx,vg_indxe).axpy(G_body,rho_g*alpha_g_old).rs();
}

}

  A_grad_U.prod(Ajac,grad_U,-1,1,-2,-1,-2);

/*
  if ((options & COMP_SOURCE) || (options & COMP_SOURCE_LUMPED)){

  if (comp_interphase_terms==1) {
    FastMat2::branch();
    // Phase interaction forces
    if (Rey_bubble>0){
        FastMat2::choose(0);

	tmp3_drag = tmp1_drag*f_drag_alpha;

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

  }  // ending COMP_SOURCE
*/

    // DEBUG

  Cp.scale(mask_Cp);
  Ajac.scale(mask_Ajac);
  Djac.scale(mask_Djac);
  Cjac.scale(mask_Cjac);

    if (0 && comp_liq_prof){
	int kk,ielhh;
	element.position(kk,ielhh);
	if (kk<10) {
	printf("Element %d \n",kk);
//        U.print("Estado :");
//    Cp.print("Cp:");
//    flux.print("flux:");
//    Ajac.print("Ajac:");
//    fluxd.print("fluxd:");
//    Djac.print("Djac:");
//    Cjac.print("Cjac:");
//    G_source.print("G_source:");
    }
    // END DEBUG
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

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  double rho_m,tau;

    const FastMat2 &grad_N = *new_adv_dif_elemset->grad_N();
    // esto es equivalente a lo viejo
    if (0) {
    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    P_supg.prod(Ao_grad_N,tau_supg_c,1,2,-1,-1,3);
    } else if(use_pmm_method && comp_gas_prof && 0) {
    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    //    P_supg.prod(Ao_grad_N,Tau_beta,1,2,-1,-1,3);
    tau=double(tau_supg_c.get(vg_indx,vg_indx));
    P_supg.set(Ao_grad_N).scale(tau);
	} else if(flag_interfacial_pressure>0 && comp_gas_prof && upwind_gas_matrix) {

/// HACER ALGO ESPECIAL
    compute_tau_gas();

//    P_supg.prod(grad_N,tau_A,-1,1,-1,2,3);

    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    P_supg.prod(Ao_grad_N,tau_A,1,2,-1,-1,3);

    }
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
    // rho_m=arho_l+arho_g;
    rho_m=rho_l;
    tau=tau/rho_m;
    tmp7.set(grad_N).scale(tau);
    tmp7.t();
    P_supg.ir(2,1).is(3,vl_indx,vl_indxe).add(tmp7);
    P_supg.rs();
    }

    if (comp_gas_prof) {

      // epsilon perturbation function for the gas phase
      tau=double(tau_supg_c.get(vg_indx,vg_indx));

      U.is(1,vg_indx,vg_indxe);
      v_g.set(U);
      U.rs();
      tmp9.prod(v_g,grad_N,-1,-1,1).scale(tau);
//      if (use_pmm_method) tmp9.scale(1./beta_pmm);
      P_supg.ir(2,2).ir(3,2).add(tmp9).rs();
      if (use_pspg_for_gas>0) {
      rho_m=rho_g;
//      if (use_pmm_method) rho_m = rho_m*beta_pmm;
      tmp7.rs();
      tmp7.set(grad_N).scale(factor_pspg_for_gas*tau/rho_m);
      tmp7.t();
      P_supg.ir(2,2).is(3,vg_indx,vg_indxe).add(tmp7).rs();
//      if (flag_debug) printf(" beta_pmm (comp P_supg )= %15.5e \n",beta_pmm);

      // added stabilization for continuity equation
      /*
	tmp7.rs();
	tmp7.set(grad_N).scale(factor_pspg_for_gas*tau*rho_m);
	tmp7.t();
	P_supg.is(2,vg_indx,vg_indxe).ir(3,2).add(tmp7);
	P_supg.rs();
      */

 }
    }

    // added stabilization for continuity equation
    if (shocap>0) {
      // stabilization enhancement of continuity equation
      tau=double(tau_supg_c.get(1,1));
      rho_m=rho_l;
      tau=tau*rho_m;
      tmp7.rs();
      tmp7.set(grad_N).scale(tau);
      tmp7.t();
      P_supg.is(2,vl_indx,vl_indxe).ir(3,1).add(tmp7);
      P_supg.rs();
    }
    tmp7.rs();

    // debug
    // P_supg.set(1.);
    // end debug
    }
}
#endif

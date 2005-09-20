//__INSERT_LICENSE__
//$Id: bubbly_asm_mix.cpp,v 1.5.28.1 2005/09/20 00:58:34 mstorti Exp $
//
//
// <<<<<<<<<<<<<<<<<< VERSION ASM >>>>>>>>>>>>>>>>>>>>>>>>
//
//


#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly_asm_mix.h"

#define bubbly_ff bubbly_asm_mix_ff
#define bubbly bubbly_asm_mix
#define bubbly_bcconv bubbly_bcconv_asm_mix

// incluido
#include "./fm2funm.h"

/*
class MyFun : public FastMat2_funm {
public:
  // No hardening at all
  double f(double l) { return sqrt(l); }
//  double f(double l) { return l; }
} my_fun;

// fin incluido
*/

extern const char * jobinfo_fields;

extern GlobParam *GLOB_PARAM;

int bubbly::ask(const char *jobinfo,int &skip_elemset) {

  // TAKE CARE WITH DUPLICATION OF CONTRIBUTIONS TO RESIDUE !!!
  // jobinfo_field = gasliq  with name() = gasliq
  // jobinfo_field = liq     with name() = liq
  // jobinfo_field = gas     with name() = gas
  // jobinfo_field = kep     with name() = liq

   skip_elemset = 1;
   DONT_SKIP_JOBINFO(comp_prof);
   DONT_SKIP_JOBINFO(comp_res);

   /*
   if ((!strcmp(jobinfo_fields,"gasliq"))
       && (!strcmp(name(),"gasliq"))) DONT_SKIP_JOBINFO(comp_res);

   if ((!strcmp(jobinfo_fields,"liq"))
       && (!strcmp(name(),"liq"))) DONT_SKIP_JOBINFO(comp_res);

   if ((!strcmp(jobinfo_fields,"gas"))
       && (!strcmp(name(),"gas"))) DONT_SKIP_JOBINFO(comp_res);

   if ((!strcmp(jobinfo_fields,"kep"))
       && (!strcmp(name(),"liq"))) DONT_SKIP_JOBINFO(comp_res);

   //   printf(" JobInfoFields  %s  Name %s \n",jobinfo_fields,name());

   */

  return 0;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
void bubbly_ff::start_chunk(int &ret_options) {
  int ierr;
  FastMat2 tmp5;

  new_adv_dif_elemset = dynamic_cast<const NewAdvDif *>(elemset);
  elemset->elem_params(nel,ndof,nelprops);
  EGETOPTDEF_ND(elemset,int,ndim,0);
  // numero de fases presentes
  EGETOPTDEF_ND(elemset,int,nphases,1);

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

  //o Direction of gravity
  EGETOPTDEF_ND(elemset,int,g_dir,ndim);

  //o stabilization factor (tau_fac)
  EGETOPTDEF_ND(elemset,double,tau_fac,1.);

  //o shocap
  EGETOPTDEF_ND(elemset,double,shocap,0.);

  //o Adjust the stability parameters, taking into account
  // the time step. If the \verb+steady+ option is in effect,
  // (which is equivalent to $\Dt=\infty$) then
  // \verb+temporal_stability_factor+ is set to 0.
  EGETOPTDEF_ND(elemset,double,temporal_stability_factor,1.);

  // bubble diameter
  EGETOPTDEF_ND(elemset,double,d_bubble,0);

  // interphase drag force term key
  EGETOPTDEF_ND(elemset,int,comp_interphase_terms,0);

  // flag masking for matrices (Cp,Ajac,Cjac,Djac)
  // 1111 todas prendidas
  EGETOPTDEF_ND(elemset,int,flag_debug,0);
  EGETOPTDEF_ND(elemset,int,mask_matrix,1111);

  // stronger coupling between continuous and disperse phases
  EGETOPTDEF_ND(elemset,int,coupled,1);

  // key to add stabilization for gas continuity equation as for incompressible
  // as suggested by a Finland researcher
  EGETOPTDEF_ND(elemset,int,use_pspg_for_gas,0);
  EGETOPTDEF_ND(elemset,double,factor_pspg_for_gas,1.0);

  // key for using alpha_grad_p as standard or using grad_p
  EGETOPTDEF_ND(elemset,int,use_alpha_grad_p,1);

  // key for using PMM method
  EGETOPTDEF_ND(elemset,int,use_pmm_method,0);

  // sonic velocity for remove singularity in Cp matrix
  EGETOPTDEF_ND(elemset,double,cc,1.e+6);

  // to activate or not the grad alpha_g in source term
  EGETOPTDEF_ND(elemset,double,flag_grad_alpha_source,1.0);

  // to scale the temporal term in the liquid continuity equation
  EGETOPTDEF_ND(elemset,double,factor_liq_mass_eq_mod,1.0);

  // to consider or not rho in disperse transport eqs. 
  EGETOPTDEF_ND(elemset,int,disperse_eqs_without_rho,0);

  // Schmidt number
  EGETOPTDEF_ND(elemset,double,Sc_number,1.0);

  //o use_modified_slag_vslip : key to modify slag slip velocity 
  //				as a function of slag void fraction 
  EGETOPTDEF_ND(elemset,int,use_modified_slag_vslip,0);

  // key for using different drag model
  EGETOPTDEF_ND(elemset,int,drag_model,0);

  //o _T: double[ndim] _N: G_body _D: null vector
  // _DOC: Vector of gravity acceleration (must be constant). _END
  G_body.resize(1,ndim);
  G_body.set(0.);
  ierr = elemset->get_double("G_body",
			     *G_body.storage_begin(),1,ndim);

  // Source for disperse phase
  EGETOPTDEF_ND(elemset,double,alpha_source,0.);
  // Fixed slip velocity imposed by user
  EGETOPTDEF_ND(elemset,double,vslip_user,0.);

  assert(ndim>0);
  assert(nphases>0);
  assert(ndof==1+ndim+nphases);

  assert(rho_l>0.);
  assert(visco_l>0.);
  assert(rho_g>0.);
  assert(visco_g>0.);

  //  phases density vector
  rho_g_vp.resize(1,nphases);
  rho_g_vp.set(rho_g);
  ierr = elemset->get_double("rho_phases",*rho_g_vp.storage_begin(),1,nphases);

  //  phases viscosity vector
  visco_g_vp.resize(1,nphases);
  visco_g_vp.set(visco_g);
  ierr = elemset->get_double("visco_phases",*visco_g_vp.storage_begin(),1,nphases);

  //  Bubble diameter vector
  d_bubble_vp.resize(1,nphases);
  d_bubble_vp.set(d_bubble);
  ierr = elemset->get_double("d_bubble_phases",*d_bubble_vp.storage_begin(),1,nphases);

  //  slip velocity vector
  vslip_user_vp.resize(1,nphases);
  vslip_user_vp.set(vslip_user);
  ierr = elemset->get_double("vslip_user_phases",*vslip_user_vp.storage_begin(),1,nphases);

  //  source term vector
  alpha_source_vp.resize(1,nphases);
  alpha_source_vp.set(alpha_source);
  ierr = elemset->get_double("alpha_source_phases",*alpha_source_vp.storage_begin(),1,nphases);

  vl_indx = 2;
  vl_indxe = vl_indx+ndim-1;
  alpha_indx = vl_indxe+1;
  alpha_indx_vp.resize(nphases);
  for (int j=0; j<nphases; j++) alpha_indx_vp[j] = vl_indxe+j+1;

  alpha_g_vp.resize(1,nphases);
  alpha_g_vp_ctf.resize(1,nphases);

  arho_g_vp.resize(1,nphases);

  v_l.resize(1,ndim);
  v_mix.resize(1,ndim);
  v_g.resize(1,ndim);
  v_g_vp.resize(2,ndim,nphases);

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

  v_rel.resize(1,ndim);
  v_rel_vp.resize(2,ndim,nphases);

  id_liq = 1, id_gas = -1;

  v_l_old.resize(1,ndim);
  v_g_old.resize(1,ndim);
  v_g_vp_old.resize(2,ndim,nphases);

  Id_vp.resize(2,nphases,nphases);
  ones_vp.resize(1,nphases);
  vslip_vp.resize(1,nphases);

  grad_alpha_g_vp.resize(2,ndim,nphases);

  visco_g_eff_vp.resize(1,nphases);
  delta_supg_vp.resize(1,nphases);

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
  alpha_g = U.get(alpha_indx);

  int istart = alpha_indx_vp[0];
  int iend   = alpha_indx_vp[nphases-1];

  U.is(1,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
  alpha_g_vp.set(U);
  U.rs();

  double alpha_d_sum = (double) alpha_g_vp.sum_all();
  alpha_l = 1.0 - alpha_d_sum;


  // limito las fracciones de vacio entre 0 y 1
  double tol=1e-5;
#if 0
  //  double alpha_g_ctf, alpha_l_ctf;

  for (int j=1; j<=nphases; j++) {
	 alpha_g = alpha_g_vp.get(j);
	 alpha_g_ctf = (alpha_g < tol ? tol : alpha_g);
	 alpha_g_vp_ctf.setel(alpha_g_ctf,j);
  }

  alpha_d_sum = (double) alpha_g_vp_ctf.sum_all();
  alpha_l = 1.0 - alpha_d_sum;

  FastMat2::branch();
  if (alpha_l<0) {
    FastMat2::choose(0);

    alpha_l_ctf = tol;
    alpha_l =alpha_l_ctf;
    alpha_g_vp_ctf.scale((1-tol)/alpha_d_sum);
  } else {
    FastMat2::choose(1);
    alpha_d_sum = 0.0;
    for (int j=1; j<=nphases; j++) {
      alpha_g_ctf = alpha_g_vp_ctf.get(j);
     alpha_g_ctf = (alpha_g_ctf > (1.0-alpha_d_sum-tol) ? (1.0-alpha_d_sum-tol) : alpha_g_ctf);
     alpha_g_vp_ctf.setel(alpha_g_ctf,j);
     alpha_d_sum += alpha_g_ctf;
    }
    alpha_l_ctf = 1.0-alpha_d_sum;
  }

  FastMat2::leave();

//  alpha_l   = alpha_l_ctf;
  //  alpha_g   = alpha_g_ctf;

  alpha_g_vp.set(alpha_g_vp_ctf);

#endif

#if 0
  double alpha_d_sum_max = 1.1;
  FastMat2::branch();
  if (alpha_d_sum>alpha_d_sum_max) {
    FastMat2::choose(0);
    alpha_g_vp.scale(1.1/alpha_d_sum);
    alpha_l = 1 - alpha_d_sum_max;
    alpha_d_sum = alpha_d_sum_max;
  }
  FastMat2::leave();
#endif

  arho_l = alpha_l*rho_l;
  arho_g = alpha_g*rho_g;

  // rho_g diagonal matrix
  Id_vp.set(0.).d(1,2);
  Id_vp.set(rho_g_vp).rs();

  arho_g_vp.prod(Id_vp,alpha_g_vp,1,-1,-1);

  // v_l is for mixture velocity
  U.is(1,vl_indx,vl_indxe);
  v_l.set(U);
  U.rs();

  double rb = d_bubble/2.;
  if (!vslip_user) {
    vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	     rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));

    vslip = (g_dir > 0 ? vslip : -vslip);
  } else vslip = vslip_user;

  if (vslip_user_vp.sum_abs_all()>0) {
	  vslip_vp.set(vslip_user_vp);
  } else {
    for (int j=1; j<=nphases; j++) {
		rb = d_bubble_vp.get(j);
        vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	     rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));

        vslip = (g_dir > 0 ? vslip : -vslip);
        vslip_vp.setel(vslip,j);
	}
  }

  // mixture density
  rho_m = arho_l+arho_g;
  rho_m = arho_l+arho_g_vp.sum_all();

  // velocidad slip respecto a la de la mezcla = a la del usuario
  vslip_m = vslip;
  vslip_m_vp.set(vslip_vp);

  // modifico velocidad slip de la escoria por la diferencia de densidades con la mezcla
  if(use_modified_slag_vslip) {
    rho_g = rho_g_vp.get(nphases);
    vslip_m_vp.ir(1,nphases).scale(1.-rho_g/rho_m).rs();
  }
  
  // velocidad del gas
  v_g.set(v_l).addel(vslip_m,abs(g_dir));

  for (int j=1; j<=nphases; j++) {
    v_g_vp.ir(2,j);
    v_g_vp.set(v_l).addel(vslip_m_vp.get(j),abs(g_dir));
    v_g_vp.rs();
  }

  FastMat2::branch();
  if (rho_m<tol) {
    FastMat2::choose(0);
    alpha_l_ctf = alpha_l;
  }
  FastMat2::leave();

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
  // Mixture
  H.setel(rho_m,1);
  H.is(1,vl_indx,vl_indxe).set(v_l).scale(rho_m);
  H.rs();
  }
  if (comp_gas_prof) {
  // Gas phase

    /*
  H.setel(arho_g,alpha_indx);
  H.rs();
    */

  H.is(1,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
  if(disperse_eqs_without_rho==0) {
    H.set(arho_g_vp);
  } else {
    H.set(alpha_g_vp);
  } 
  H.rs();

  }

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
void bubbly_ff::compute_delta_sc_v(FastMat2 &delta_sc_v) {

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  // Improve it latter, please ....

  delta_sc_v.set(0.);
  if (comp_gas_prof) {
    delta_sc_v.is(1,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
    //    delta_sc_v.set(delta_supg).rs();
    delta_sc_v.set(delta_supg_vp).rs();
  }

}


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

    double Peclet = velmod * h_supg / (2. * visco_supg);
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

  int comp_liq_prof  = (!strcmp(jobinfo_fields,"liq")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_gas_prof  = (!strcmp(jobinfo_fields,"gas")|!strcmp(jobinfo_fields,"gasliq"));
  int comp_kep_prof  = !strcmp(jobinfo_fields,"kep");

  //  Matrix seed;
  //  seed= Matrix(ndof,ndof);
  seed.set(0.);

  if (comp_liq_prof && comp_gas_prof) {
    for (int j=1; j<=ndim+1+nphases; j++) {
      for (int k=1; k<=ndim+1+nphases; k++) {
        seed.setel(1.,j,k);
      }
    }
  }
  else if (comp_liq_prof) {
//    int ip[] = {1,2,3,4};
    for (int j=1; j<=ndim+1; j++) {
      for (int k=1; k<=ndim+1; k++) {
        seed.setel(1.,j,k);
      }
    }
  }

  else if (comp_gas_prof) {
//    int ip[] = {2};
    for (int j=ndim+2; j<=ndim+1+nphases; j++) {
      for (int k=ndim+2; k<=ndim+1+nphases; k++) {
        seed.setel(1.,j,k);
      }
    }
  }

  /*
  int nelprops,nel,ndof;
  elemset->elem_params(nel,ndof,nelprops);

  vector<double> bcconv_factor;
  const char *line;
  elemset->get_entry("bcconv_factor",line);
  if(line) {
    read_double_array(bcconv_factor,line);
    PETSCFEM_ASSERT0(bcconv_factor.size() == ndof,
		     "bcconv_factor needs ndof values \n");
  } else bcconv_factor.resize(ndof,1.0);

  for (int j=0; j<ndof; j++) {
	 seed.ir(1,j+1).scale(bcconv_factor[j]).rs();
	 seed.ir(2,j+1).scale(bcconv_factor[j]).rs();
	  }
  */

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

  //  int alpha_indx = ndim+1+1;
  Cp.set(0.);
  if (comp_liq_prof) {

  // mixture momentum conservation
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(rho_m).rs();

  if(coupled) {
    for (int j=1; j<=nphases; j++) {
      rho_g = rho_g_vp.get(j);
  // mixtures mass conservation
  Cp.setel(rho_g-rho_l,1,alpha_indx_vp[j-1]);
  // mixture momentum conservation
  Cp.is(1,vl_indx,vl_indxe).ir(2,alpha_indx_vp[j-1]).set(v_l).scale(rho_g-rho_l).rs();
  }
  }
  }

  if (comp_gas_prof) {
  // gas mass conservation
    /*
  Cp.setel(rho_g,alpha_indx,alpha_indx);
    */
  Id_vp.set(0.).d(1,2);
  
  if (disperse_eqs_without_rho==0) {
  Id_vp.set(rho_g_vp).rs();
  } else {
  Id_vp.set(1.0).rs();
  }

  Cp.is(1,alpha_indx_vp[0],alpha_indx_vp[nphases-1]).is(2,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
  Cp.set(Id_vp).rs();

  }

  Cp.scale(mask_Cp);
  // Cp.print("Cp:\n");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Advective fluxes
  Amoml.prod(v_l,v_l,1,2).scale(rho_m).axpy(Id,p);
  double vaux; 
  
  flux.set(0.);
  if (comp_liq_prof) {
    // mixture mass conservation
    flux.ir(1,1).set(v_l).scale(rho_m).rs();
    // mixture momentum conservation
    flux.is(1,vl_indx,vl_indxe).set(Amoml).rs();

    flux.ir(1,vl_indx-1+abs(g_dir)).ir(2,abs(g_dir));
    for (int j=1; j<=nphases; j++) {
      rho_g = rho_g_vp.get(j);
      alpha_g = alpha_g_vp.get(j);
      vslip_m = vslip_m_vp.get(j);
      vaux = rho_g*alpha_g*vslip_m*vslip_m;
      flux.add(vaux);
    }
    flux.rs();
  }
  
  if (comp_gas_prof) {
    // gas mass conservation
    /*
  flux.ir(1,alpha_indx).set(v_g).scale(arho_g).rs();
    */

  Id_vp.set(0.).d(1,2);
  if (disperse_eqs_without_rho==0) {
  Id_vp.set(arho_g_vp).rs();
  } else {
  Id_vp.set(alpha_g_vp).rs();
  }

  flux.is(1,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
  flux.prod(Id_vp,v_g_vp,1,-1,2,-1).rs();

  }
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Adjective Jacobians
  Ajac.set(0.);
  
  if (comp_liq_prof) {
    // mixture mass conservation
    Ajac.ir(2,1).is(3,vl_indx,vl_indxe).axpy(Id,rho_m).rs();
    
    if(coupled) {
      for (int j=1; j<=nphases; j++) {
	rho_g = rho_g_vp.get(j);
	Ajac.ir(2,1).ir(3,alpha_indx_vp[j-1]).set(v_l).scale(rho_g-rho_l).rs();
	Ajac.is(2,vl_indx,vl_indxe).ir(3,alpha_indx_vp[j-1]).prod(v_l,v_l,1,2).scale(rho_g-rho_l).rs();
	vslip_m = vslip_m_vp.get(j);
	vaux = rho_g*vslip_m*vslip_m;
	Ajac.ir(1,abs(g_dir)).is(2,vl_indx-1+abs(g_dir)).ir(3,alpha_indx_vp[j-1]).add(vaux).rs();
      }
    }
    
    // mixture momentum conservation
    Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).axpy(Id,1.0);
    
    Y.rs().prod(v_l,Id,2,1,3).scale(rho_m);
    Ajac.rs().is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe).add(Y);
    Y.prod(v_l,Id,1,2,3).scale(rho_m).rs();
    Ajac.add(Y);
    Y.rs();
    Ajac.rs();
  }
  
  if (comp_gas_prof) {
    // gas mass conservation
    
    for (int j=1; j<=nphases; j++) {
      v_g_vp.ir(2,j);
      rho_g = rho_g_vp.get(j);
      alpha_g = alpha_g_vp.get(j);
      // alpha_l = 1-alpha_g;
      Ajac.ir(2,alpha_indx_vp[j-1]).ir(3,alpha_indx_vp[j-1]);
      Ajac.set(v_g_vp);
      if (disperse_eqs_without_rho==0) Ajac.scale(rho_g);

      v_g_vp.rs();
      Ajac.rs();
      
      /*
	vslip = vslip_vp.get(j);
	v_rel.set(v_l).scale(0.0);
	v_rel.addel(vslip,abs(g_dir));
      */
      if(coupled) {	
	arho_g = arho_g_vp.get(j);
	if (disperse_eqs_without_rho==0) {
	Ajac.ir(2,alpha_indx_vp[j-1]).is(3,vl_indx,vl_indxe).axpy(Id,arho_g).rs();
	} else {
	Ajac.ir(2,alpha_indx_vp[j-1]).is(3,vl_indx,vl_indxe).axpy(Id,alpha_g).rs();
	}	
      }
    }    
    /*  
	if(coupled && use_modified_slag_vslip) {	
	// agregado termino del jacobiano por modificacion de la velocidad slip en la escoria
	rho_g = rho_g_vp.get(nphases);
	arho_g = arho_g_vp.get(nphases);
	alpha_g = alpha_g_vp.get(nphases);
	vslip = vslip_vp.get(nphases);
	if (disperse_eqs_without_rho==0) {
	vaux = arho_g*vslip*(rho_g-rho_l)*rho_g/rho_m/rho_m;
	Ajac.ir(1,abs(g_dir)).ir(2,alpha_indx_vp[nphases-1]).ir(3,alpha_indx_vp[nphases-1]).add(vaux).rs();
	} else {
	vaux = alpha_g*vslip*(rho_g-rho_l)*rho_g/rho_m/rho_m;
	Ajac.ir(1,abs(g_dir)).ir(2,alpha_indx_vp[nphases-1]).ir(3,alpha_indx_vp[nphases-1]).add(vaux).rs();
	}
	}
    */
  }
  Ajac.scale(mask_Ajac);
  
  // mixture mass conservation
  // mixture momentum conservation
  // gas mass conservation

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

//  grad_U.ir(2,alpha_indx);
//  grad_alpha_g.set(grad_U);
//  grad_U.rs();
//  grad_alpha_l.set(grad_alpha_g).scale(-1.0);

  // grad_alpha for phases vector
  grad_U.is(2,alpha_indx_vp[0],alpha_indx_vp[nphases-1]);
  grad_alpha_g_vp.set(grad_U);
  grad_U.rs();
  grad_alpha_l_vp.set(grad_alpha_g_vp).scale(-1.0);

  grad_U.ir(2,1);
  grad_p.set(grad_U);
  grad_U.rs();

  // Smagorinsky turbulence model
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

  visco_m_eff = alpha_l*visco_l+alpha_g*visco_g + visco_sato + rho_m*visco_t;
  visco_g_eff = rho_g*visco_t/Sc_number;

  visco_m_eff = (visco_m_eff <= 0 ? 1.0e-15 : visco_m_eff);
  visco_g_eff = (visco_g_eff <= 0 ? 1.0e-15 : visco_g_eff);

  fluxd.set(0.);
  Djac.set(0.);

  if (comp_liq_prof) {
    // mixture phase
    fluxd.is(1,vl_indx,vl_indxe).set(strain_rate_l).scale(2.*visco_m_eff).rs();
    Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,visco_m_eff).rs();
  }

  if (comp_gas_prof) {
    // Gas phase

    /*
      fluxd.ir(1,alpha_indx).set(grad_alpha_g).scale(visco_g_eff).rs();
      Djac.ir(2,alpha_indx).ir(4,alpha_indx).axpy(Id,visco_g_eff).rs();
    */

    for (int j=1; j<=nphases; j++) {
      rho_g = rho_g_vp.get(j);
      grad_alpha_g_vp.ir(2,j);
      grad_alpha_g.set(grad_alpha_g_vp);
      visco_g_eff = visco_t/Sc_number;
      if (disperse_eqs_without_rho==0) visco_g_eff = visco_g_eff * rho_g;
      visco_g_eff = (visco_g_eff <= 0 ? 1.0e-15 : visco_g_eff);
      visco_g_eff_vp.setel(visco_g_eff,j);
      fluxd.ir(1,alpha_indx_vp[j-1]).set(grad_alpha_g).scale(visco_g_eff).rs();
      Djac.ir(2,alpha_indx_vp[j-1]).ir(4,alpha_indx_vp[j-1]).axpy(Id,visco_g_eff).rs();
      grad_alpha_g_vp.rs();
    }
  }

  Djac.scale(mask_Djac);

  G_source.set(0.);
  Cjac.set(0.);
  
  if (options & COMP_SOURCE) {
    // gravity forces
    if (comp_liq_prof) {
      G_source.is(1,vl_indx,vl_indxe).axpy(G_body,rho_m).rs();
      
      if(coupled) {	
	for (int j=1; j<=nphases; j++) {
	  rho_g = rho_g_vp.get(j);
	  Cjac.is(1,vl_indx,vl_indxe).ir(2,alpha_indx_vp[j-1]).axpy(G_body,-(rho_g-rho_l)).rs();
	}
      }
    }
    
    
    if (comp_gas_prof) {
      
      //    G_source.addel(alpha_source*rho_g,alpha_indx);
      for (int j=1; j<=nphases; j++) {
	rho_g = rho_g_vp.get(j);
	alpha_source = alpha_source_vp.get(j);
	if (disperse_eqs_without_rho==0) {
	G_source.addel(alpha_source*rho_g,alpha_indx_vp[j-1]);
	} else {
	G_source.addel(alpha_source,alpha_indx_vp[j-1]);
	}
      }
    }
  }

  Cjac.scale(mask_Cjac);

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
    v_g_old.set(v_l_old).addel(vslip_m,abs(g_dir));

    for (int j=1; j<=nphases; j++) {
      v_g_vp_old.ir(2,j);
      v_g_vp_old.set(v_l_old).addel(vslip_m_vp.get(j),abs(g_dir));
      v_g_vp_old.rs();
    }

    const FastMat2 &grad_N = *advdf_e->grad_N();
    vel_supg.set(v_l_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }

    visco_supg = visco_m_eff/rho_m;
    int ijob =0;
    compute_tau(ijob);

    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }
    tau_supg.eye(tau_supg_a).setel(delta_supg,1,1);

    delta_supg_vp.set(0.);

    for (int j=1; j<=nphases; j++) {

      v_g_old.set(v_g_vp_old.ir(2,j));
      rho_g = rho_g_vp.get(j);

      vel_supg.set(v_g_old);
      if(axi>0){
	vel_supg.setel(0.,axi);
      }
      
      visco_g_eff = visco_g_eff_vp.get(j);
      if (disperse_eqs_without_rho==0) { 
	visco_supg = visco_g_eff/rho_g;
      } else {
	visco_supg  = visco_g_eff;
      }
      
      ijob=0;
      compute_tau(ijob);

      if (tau_fac != 1.) {
	tau_pspg *= tau_fac;
	tau_supg_a *= tau_fac;
      }

      tau_supg.setel(tau_supg_a,alpha_indx_vp[j-1],alpha_indx_vp[j-1]);
      /*      
      if (disperse_eqs_without_rho==0) {
	//      delta_supg_vp.setel(delta_supg,j).scale(rho_g);
      delta_supg_vp.ir(1,j).set(delta_supg).scale(rho_g).rs();

      } else {
	delta_supg_vp.setel(delta_supg,j);
      //      printf(" delta_supg : phase # %d valor: %f \n",j,delta_supg);
      */
      delta_supg_vp.setel(delta_supg,j);
      
      v_g_vp_old.rs();

    }
    tau_supg_c.set(tau_supg);
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

  double tau;

  const FastMat2 &grad_N = *new_adv_dif_elemset->grad_N();
  // esto es equivalente a lo viejo

  //    rho_m = alpha_g*rho_g+alpha_l*rho_l;

  if (0) {
    const FastMat2 &Ao_grad_N = new_adv_dif_elemset->Ao_grad_N;
    P_supg.prod(Ao_grad_N,tau_supg_c,1,2,-1,-1,3);
  } else {
    P_supg.set(0.);
    U.is(1,vl_indx,vl_indxe);
    v_l.set(U);
    U.rs();
    if (comp_liq_prof) {

      // delta perturbation function for the liquid phase
      tau=double(tau_supg_c.get(vl_indx,vl_indx));
      tmp9.prod(v_l,grad_N,-1,-1,1).scale(tau);
      tmp6.prod(tmp9,Id,1,2,3);

      P_supg.is(2,vl_indx,vl_indxe).is(3,vl_indx,vl_indxe);
      P_supg.add(tmp6).rs();

      // epsilon perturbation function for the liquid phase
      //      tau=tau/rho_m;
      tau=tau/rho_l;

      tmp7.set(grad_N).scale(tau);
      tmp7.t();
      P_supg.ir(2,1).is(3,vl_indx,vl_indxe).add(tmp7).rs();

      // added stabilization for continuity equation
      if (shocap>0 ) {
	// stabilization enhancement of continuity equation
	tau=double(tau_supg_c.get(1,1));
	//	tau=tau*rho_m;
	tau=tau*rho_l;
	tmp7.rs();
	tmp7.set(grad_N).scale(tau);
	tmp7.t();
	P_supg.is(2,vl_indx,vl_indxe).ir(3,1).add(tmp7);
	P_supg.rs();
      }
      tmp7.rs();

    }

    if (comp_gas_prof) {
      // supg perturbation function for the gas phase
      for (int j=1; j<=nphases; j++) {
	tau=double(tau_supg_c.get(alpha_indx_vp[j-1],alpha_indx_vp[j-1]));
	//      v_g.set(v_l).addel(vslip_m,abs(g_dir));
	v_g.set(v_g_vp.ir(2,j));
	v_g_vp.rs();
	tmp9.prod(v_g,grad_N,-1,-1,1).scale(tau);
	P_supg.ir(2,alpha_indx_vp[j-1]).ir(3,alpha_indx_vp[j-1]).add(tmp9).rs();
      }
    }

    // debug
    // P_supg.set(1.);
    // end debug
  }
}
#endif

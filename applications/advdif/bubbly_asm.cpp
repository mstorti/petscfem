//__INSERT_LICENSE__
//$Id: bubbly_asm.cpp,v 1.1 2004/09/30 16:52:35 mstorti Exp $
//
//
// <<<<<<<<<<<<<<<<<< VERSION ASM >>>>>>>>>>>>>>>>>>>>>>>>
//
//

#include <src/fem.h>
#include <src/texthash.h>
#include <src/getprop.h>

#include "bubbly_asm.h"

#define bubbly_ff bubbly_asm_ff
#define bubbly bubbly_asm
#define bubbly_bcconv bubbly_bcconv_asm

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

  // flag masking for matrices (Cp,Ajac,Cjac,Djac)
  // 1111 todas prendidas
  EGETOPTDEF_ND(elemset,int,flag_debug,0);
  EGETOPTDEF_ND(elemset,int,mask_matrix,1111);

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

  // Schmidt number
  EGETOPTDEF_ND(elemset,double,Sc_number,1.0);

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
  assert(ndof==2+ndim);
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
	  rho_l_eff= rho_l*(1.0-C_vm);
	  rho_g_eff= rho_g + C_vm*rho_l;
	  tmp_vm.resize(2,ndim,ndim);
	  tmp2_vm.resize(1,ndim);
	  Dmat_vl.resize(1,ndim);
	  Dmat_vg.resize(1,ndim);
  } else {
	  rho_l_eff = rho_l;          ;
	  rho_g_eff = rho_g;
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

#if 0
  // limito las fracciones de vacio entre 0 y 1
  double tol=1e-4;
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

  // Velocities
  U.is(1,vl_indx,vl_indxe);
  v_l.set(U);
  U.rs();

  if (!vslip_user) {
    double rb = d_bubble/2.;
    vslip = (rb<7e-4 ? 4474*pow(rb,1.357) :
	     rb<5.1e-3 ? 0.23 : 4.202*pow(rb,0.547));
    
    vslip = (g_dir > 0 ? vslip : -vslip);
  } else vslip = vslip_user;
  v_g.set(v_l).addel(vslip,abs(g_dir));
  
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
  if (comp_gas_prof) {
  // Gas phase
  H.setel(arho_g,2);
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

  delta_sc_v.set(0.);
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
    for (int j=1; j<=ndim+2; j++) {
      for (int k=1; k<=ndim+2; k++) {
        seed.setel(1.,j,k);
      }
    }
  }
  else if (comp_liq_prof) {
    int ip[] = {1,3,4,5};
    for (int j=0; j<=ndim; j++) {
      for (int k=0; k<=ndim; k++) {
        seed.setel(1.,ip[j],ip[k]);
      }
    }
  }

  else if (comp_gas_prof) {
    int ip[] = {2};
    for (int j=0; j<1; j++) {
      for (int k=0; k<1; k++) {
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
  // Cp.is(1,vl_indx,vl_indxe).ir(2,2).set(v_l).scale(-rho_l_eff).rs();
  Cp.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).eye(arho_l_eff).rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  Cp.setel(rho_g,2,2);
  }
  // Cp.print("Cp:\n");

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Advective fluxes
  flux.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  flux.ir(1,1).set(v_l).scale(arho_l).rs();
  // start MAR 02 2004 ==>
  // Amoml.prod(v_l,v_l,1,2).scale(rho_l_eff).axpy(Id,p);
  Amoml.prod(v_l,v_l,1,2).scale(rho_l_eff);
  // end MAR 02 2004 ==>

  if (use_alpha_grad_p) {
    Amoml.axpy(Id,p);
    flux.is(1,vl_indx,vl_indxe).axpy(Amoml,alpha_l).rs();
  } else {
    flux.is(1,vl_indx,vl_indxe).axpy(Amoml,alpha_l).rs();
    flux.is(1,vl_indx,vl_indxe).axpy(Id,p).rs();
  }
  }
  
  if (comp_gas_prof) {
  // Gas phase
  flux.ir(1,2).set(v_g).scale(arho_g).rs();
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  // Adjective Jacobians
  Ajac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Ajac.rs().ir(2,1).ir(3,2).set(v_l).scale(-rho_l);
  Ajac.rs().ir(2,1).is(3,vl_indx,vl_indxe).axpy(Id,arho_l);
//  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,2).axpy(Amoml,-1);

 // por sugerencia de Carlos introducimos una modificacion ya que deberiamos
 // estar resolviendo la ecuacion de momento de la mezcla entonces eld alpha_l
 // lo reemplazamos por la unidad
  if (use_alpha_grad_p) {
  Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).axpy(Id,alpha_l);
} else {
   Ajac.rs().is(2,vl_indx,vl_indxe).ir(3,1).axpy(Id,1.0);
}

  if (0) {
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
  //
  // Gas phase
  // aqui puede estar faltando el termino que tiene cuenta que la ec continuidad del
  // gas depende de la velocidad del liquido ya que la vel_gas=vel_liq+vel_slip
  //
  Ajac.ir(2,2).ir(3,2).set(v_g).scale(rho_g).rs();
  Ajac.ir(2,2).is(3,vl_indx,vl_indxe).axpy(Id,arho_g);
  Ajac.rs();
  }

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

  grad_U.ir(2,2);
  grad_alpha_g.set(grad_U);
  grad_U.rs();
  grad_alpha_l.set(grad_alpha_g).scale(-1.0);

  grad_U.ir(2,1);
  grad_p.set(grad_U);
  grad_U.rs();

if (0) {
  // Turbulent viscosity
  // fixme:= la tasa de deformacion se toma la del liquido??
  strain_rate_scalar = strain_rate_l.sum_square_all();
  // fixme:= esto no compila
  visco_t = C_mu*rho_l * square(k)/eps;
  P_k = 2*visco_t*strain_rate_scalar;
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
  visco_g_eff = visco_g_eff/Sc_number;

  fluxd.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  fluxd.is(1,vl_indx,vl_indxe).set(strain_rate_l).scale(2.*visco_l_eff).rs();
  }

  if (comp_gas_prof) {
  // Gas phase
  fluxd.ir(1,2).set(grad_alpha_g).scale(visco_g_eff).rs();
  }


  Djac.set(0.);
  if (comp_liq_prof) {
  // Liquid phase
  Djac.is(2,vl_indx,vl_indxe).is(4,vl_indx,vl_indxe).axpy(IdId,visco_l_eff).rs();
  }

  if (comp_gas_prof) {
  // Liquid phase
  Djac.ir(2,2).ir(4,2).axpy(Id,visco_g_eff).rs();
  }

  // Reactive terms
  //   [1] bouyancy terms
  //   [2] Phase interaction forces
  //   [3] pressure gradient terms for gas phase

  G_source.set(0.);
  Cjac.set(0.);

  if ((options & COMP_SOURCE) || (options & COMP_SOURCE_NOLUMPED)){

  // gravity forces
  if (comp_liq_prof) {
    Cjac.is(1,vl_indx,vl_indxe).ir(2,2).axpy(G_body,rho_l).rs();
    G_source.is(1,vl_indx,vl_indxe).axpy(G_body,arho_l).rs();
  }
  if (comp_gas_prof) {
    G_source.addel(alpha_source*rho_g,2);
    // static int kk=0;
    // if (!(kk++ % 500)) printf("elemset %p, agregando alpha_source %f\n",
    // elemset,alpha_source);
  }


    if (comp_liq_prof) {
    advdf_e = dynamic_cast<const NewAdvDif *>(elemset);
    assert(advdf_e);

    FastMat2 &Uo = (FastMat2 &) advdf_e->Uold();
    FastMat2 &grad_Uo = (FastMat2 &) advdf_e->grad_Uold();

    grad_Uo.ir(2,2);
    grad_alpha_g_old.set(grad_Uo);
    grad_Uo.rs();
    grad_alpha_l_old.set(grad_alpha_g_old).scale(-1.0);
    double pold = Uo.get(1);

    // G_source.is(1,vl_indx,vl_indxe).axpy(grad_alpha_l_old,pold*flag_grad_alpha_source).rs();
    // G_source.is(1,vl_indx,vl_indxe).axpy(grad_p,-alpha_l*flag_grad_alpha_source).rs();

  // virtual mass acceleration for liquid phase
    if (comp_virtual_mass>0) {

      Uo.is(1,vl_indx,vl_indxe);
	  v_l_old.set(Uo);
	  Uo.rs();

	  v_g_old.set(v_l_old).addel(vslip,abs(g_dir));

	  Dmat_vl.prod(v_l,grad_v_l,-1,1,-1).axpy(v_l,rec_Dt).axpy(v_l_old,-rec_Dt);
	  Dmat_vg.prod(v_g,grad_v_g,-1,1,-1).axpy(v_g,rec_Dt).axpy(v_g_old,-rec_Dt);

	  tmp2_vm.set(Dmat_vg).scale(alpha_g).rest(Dmat_vl);
	  tmp2_vm.scale(rho_l*C_vm);
	  G_source.is(1,vl_indx,vl_indxe).add(tmp2_vm).rs();

	  tmp_vm.set(grad_v_l).axpy(Id,rec_Dt).scale(-rho_l*C_vm);
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).add(tmp_vm).rs();

}
}
}

  if ((options & COMP_SOURCE) || (options & COMP_SOURCE_LUMPED)){

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
      if (drag_model==0) {
      C_drag_ff= 24.0/Rey_bubble;
      dCDdRe_ff= -24.0/Rey_bubble/Rey_bubble;
      } else if (drag_model==1) {
      double tmp0_drag = (1.0+0.15*pow(Rey_bubble,0.687));
	  C_drag_ff= 24.0/Rey_bubble*tmp0_drag;
	  dCDdRe_ff= -24.0/Rey_bubble/Rey_bubble*tmp0_drag +
                      24.0/Rey_bubble*(0.15*0.687*pow(Rey_bubble,0.687-1.0));
  } else if (drag_model==2) {
	  C_drag_ff= 0.44;
	  dCDdRe_ff= 0.0;
  }

	  dRedU = rho_l*d_bubble/visco_l;

	  tmp1_drag = C1_drag*Rey_bubble*C_drag_ff;

      coe = 1.0;
      pot = 1.0;
      f_drag_alpha = coe*pow(alpha_l,pot)*alpha_g;
      df_drag_alpha = coe*pow(alpha_l,pot-1)*(1.-(1.+pot)*alpha_g);

      tmp2_drag = tmp1_drag*df_drag_alpha;

	  tmp3_drag = tmp1_drag*f_drag_alpha;

          tmp4_drag = C_drag_ff+Rey_bubble*dCDdRe_ff;

          Phi_1.prod(v_g_l,v_g_l,1,2).scale(-tmp4_drag*dRedU/v_slip);
          Phi_2.set(Id).scale(Rey_bubble*C_drag_ff);
          Phi_1.rest(Phi_2);

	  Phi_1.scale(C1_drag*f_drag_alpha);

          if (comp_liq_prof) {
          // Liquid phase
	  Cjac.is(1,vl_indx,vl_indxe).is(2,vl_indx,vl_indxe).axpy(Phi_1,id_liq*id_liq).rs();
	G_source.is(1,vl_indx,vl_indxe).axpy(v_g_l,id_liq*tmp3_drag).rs();
          }
    }


    FastMat2::leave();
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

    v_g_old.set(v_l_old).addel(vslip,abs(g_dir));

    const FastMat2 &grad_N = *advdf_e->grad_N();

    //    double tau_supg_a,tau_pspg,delta_supg,visco_supg,velmod, h_supg;
    // Fase liquida
    vel_supg.set(v_l_old);
    if(axi>0){
      vel_supg.setel(0.,axi);
    }

    //    if (use_pmm_method) compute_tau_beta();

    visco_supg = visco_l_eff/rho_l;
    int ijob =0;
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
    ijob=0;
    compute_tau(ijob);

    if (tau_fac != 1.) {
      tau_pspg *= tau_fac;
      tau_supg_a *= tau_fac;
    }

    tau_supg.setel(tau_supg_a,2,2);

    tau_supg_c.set(tau_supg);

  }

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

#if 0
    static int kk=0;
    if (!(kk++ % 400)) {
      printf("elemset %p, G_source: \n",elemset);
      G_source.print("");
    }
#endif
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
      tau=double(tau_supg_c.get(2,2));

      v_g.set(v_l).addel(vslip,abs(g_dir));
      tmp9.prod(v_g,grad_N,-1,-1,1).scale(tau);
      P_supg.ir(2,2).ir(3,2).add(tmp9).rs();
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

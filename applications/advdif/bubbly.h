// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
// $Id: bubbly.h,v 1.14 2004/12/21 12:20:37 mstorti Exp $

#include "./advective.h"
#include "./stream.h"

#include "nwadvdifj.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
/** Flux function for multi-phase flow
 */
class bubbly_ff : public AdvDifFFWEnth {
private:
  int nelprops,nel,ndof,ndim,k_indx,e_indx,vg_indx,vl_indx,
    vl_indxe,vg_indxe;
  int coupled, disperse_eqs_without_rho;
  FastMat2 U,v_l,v_g,v_mix,Cp,Ajac,Id,Amoml,Amomg,Y,
    Djac,tmp1,Cjac,tmp2,tmp3,grad_v_l,strain_rate_l,
    grad_v_g,strain_rate_g,grad_k,grad_e,IdId,G_body,
    uintri,svec,tmp9,W_N,grad_alpha_g,grad_alpha_l,grad_p,tmp6,tmp7,tmp10;
  double alpha_source;
  FastMat2 Cpc,Ajacc,Djacc,Cjacc;
  FastMat2 Phi_1,Phi_2,v_g_l;
  FastMat2 tau_supg_c, vel_supg;
  FastMat2 tmp_debug,tmp_vm,tmp2_vm,Dmat_vl,Dmat_vg;
  FastMat2 grad_alpha_g_old,grad_alpha_l_old;
  FastMat2 grad_v_g_old,grad_v_l_old;
  FastMat2 ugn,uln;
  FastMat2 sign_Lambda_gas, VV_gas,VVi_gas;  // eigensystem decomposition of gas block
  FastMat2 v_mesh,grad_v_mesh,v_l_mesh,v_g_mesh_vp;

  int comp_virtual_mass;
  double C_vm;
  double arho_l_eff,arho_g_eff;
  double rho_l_eff,rho_g_eff;

  double alpha_l,alpha_g,arho_l,arho_g,p,k,eps,
    visco_l,visco_g,visco_l_eff,visco_g_eff,
    C_mu,C_1,C_2,sigma_k,sigma_e,P_k,tau_fac,shocap,
    visco_t,temporal_stability_factor;
  double A_van_Driest,C_smag;
  int LES, g_dir;
  FastMat2 tmp15;
  double vslip, vslip_user, Sato_model_coef, vslip_m;
  double Sc_number;

  double tmp1_drag,tmp2_drag,tmp3_drag,tmp4_drag,C1_drag,v_slip,Rey_bubble;
  double dRedU,C_drag_ff,dCDdRe_ff,id_liq,id_gas;
  double rho_l,rho_g;
  double d_bubble;
  double factor_alpha_liq;
  double flag_grad_alpha_source , flag_interfacial_pressure;
  double factor_liq_mass_eq_mod;
  int flag_predictor_scheme;

  int comp_interphase_terms,mask_matrix,use_pspg_for_gas,use_pmm_method,use_alpha_grad_p;
  int drag_model, upwind_gas_matrix;

  FastMat2 tmp1_lift,tmp2_lift,rotor_v_l;
  double drag_value , C_lift , coef_turbulent_dispersion,Sc_t;

  int flag_debug;
  double factor_pspg_for_gas, factor_sonic_speed;

  FastMat2 v_l_old,v_g_old;
  ElementIterator element;

  const NewAdvDif *advdf_e;

  //  int axi;

  void compute_tau(int ijob);

  double tau_supg_a,tau_pspg,delta_supg,visco_supg,velmod, h_supg, h_pspg;
  double cc_supg;

  void compute_beta_pmm();

  void compute_tau_beta();

  void compute_tau_gas();

  FastMat2 tmp_beta_1,tmp_beta_2,tmp_beta_3,tmp_beta_31,tmp_beta_32,tmp_beta_4,tmp_beta_40,
	  tmp_beta_41,tmp_beta_42,v_rel;
  FastMat2 Tau_beta,fTau_beta,Cpi;
  double beta_pmm, cc;

  FastMat2 temp_tau_A, temp2_tau_A, tau_A;

  // added for ASM based on mixture velocity
  double rho_m, visco_m_eff ;
  int alpha_indx , nphases;

  // extended to N phases
  vector<int> alpha_indx_vp;
  FastMat2 alpha_g_vp,alpha_g_vp_ctf,rho_g_vp,visco_g_vp,v_g_vp,v_g_vp_old,v_rel_vp;
  FastMat2 arho_g_vp,Id_vp,vslip_vp,vslip_user_vp,d_bubble_vp,vslip_m_vp,ones_vp,grad_alpha_g_vp,grad_alpha_l_vp;
  FastMat2 alpha_source_vp,visco_g_eff_vp,delta_supg_vp;
  double alpha_l_ctf, alpha_g_ctf;

public:
  bubbly_ff(NewElemset *elemset_);

  ~bubbly_ff();

  void set_profile(FastMat2 &seed);

  /** This is called before any other in a loopnike and may help in
      optimization
      @param ret_options (input/output) this is used by the flux
      function writer for returning some options. Currently the only
      option used is #SCALAR_TAU#. This options tells the elemset
      whether the flux function returns a scalar or matrix
      #tau_supg#.
  */
  void start_chunk(int &ret_options);

  /** This is called before entering the Gauss points loop and may
      help in optimization.
      @param element (input) an iterator on the elemlist.
  */
  void element_hook(ElementIterator &element);

  /** Basically stores `U(1)' in the water depth `u' and computes
      geometric parameters of the channel
      @param U (input) the state of the fluid
      @param grad_U (input) gradient of the state of the fluid
  */
  void set_state(const FastMat2 &U,const FastMat2 &grad_U);

  /** Basically stores `U(1)' in the water depth `u' and computes
      geometric parameters of the channel
      @param U (input) the state of the fluid
  */
  void set_state(const FastMat2 &U);

  /** @name Advective jacobians related */
  //@{
  /** Computes the product #(A_grad_N)_(p,mu,nu) = A_(i,mu,nu) (grad_N)_(i,p)#
      @param A_grad_N (output, size #nel# x #nd# x #nd#)
      @param grad_N (input, size #nel# x #ndof#)
  */
  void comp_A_grad_N(FastMat2 & A_grad_N,FastMat2 & grad_N);

  /** Computes the product #(A_jac_n)_(mu,nu) = A_(i,mu,nu) normal_i#
      @param A_jac_n (output, size #ndof# x #ndof#)
      @param normal (input, size #ndim#)
  */
  void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal);

  /** Computes fluxes, upwind parameters etc...
      fixme:= include more doc here ...
  */
  void compute_flux(COMPUTE_FLUX_ARGS);
  //@}

  /** @name Diffusive jacobians related */
  //@{
  /** Computes the product #(grad_N_D_grad_N)_(p,mu,q,nu)
      = D_(i,j,mu,nu) (grad_N)_(i,p) (grad_N)_(j,q)#
      @param grad_N_D_grad_N (output) size #nel# x #ndof# x #nel# x #ndof#
      @param grad_N (input) size #nel# x #ndof#
  */
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);
  //@}

  /** @name Reactive jacobians related */
  //@{
  /** Computes the product #(N_N_C)_(p,mu,q,nu)
      = w C_(mu,nu) N_p N_q#
      @param N_N_C (output, size #nel# x #ndof# x #nel# x #ndof#
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */
  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);

  /** Computes the product #(N_P_C)_(mu,q,nu)
      = w (P_supg)_(mu,lambda) C_(lambda,nu) N_q #
      @param N_P_C (output) size  #ndof# x #nel# x #ndof#
      @param P_supg (input) SUPG perturbation function size #ndof# x #ndof#
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w);

  /// Computes the enthalpy
  void enthalpy(FastMat2 &H);
  /** Computes the product of the enthalpy jacobian matrix
      with a shape function and a weight function.
  */
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  /** Computes the product of the enthalpy jacobian
      matrix with an SUPG weight function
  */
  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

#define USE_COMP_P_SUPG
#ifdef USE_COMP_P_SUPG
  void comp_P_supg(FastMat2 &P_supg);
#endif

  void get_C(FastMat2 &C);

  void get_Cp(FastMat2 &Cp);

  void compute_delta_sc_v(FastMat2 &delta_sc_v);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
/// The elemset corresponding to the `bubbly_ff' flux function.
class bubbly : public NewAdvDif {
public:
  /** Constructor, creates the fluc function object.
      fixme:= should destroy the flux functin.
  */
  bubbly() :  NewAdvDif(new bubbly_ff(this)) {};

  int ask(const char *jobinfo,int &skip_elemset);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
class bubbly_bcconv : public NewBcconv {
public:
  bubbly_bcconv() : NewBcconv(new bubbly_ff(this)) {};
};



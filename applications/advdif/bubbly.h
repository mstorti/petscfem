// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: bubbly.h,v 1.8 2002/06/26 20:59:01 mstorti Exp $
#ifndef BUBBLY_H
#define BUBBLY_H

#include "./advective.h"
#include "./stream.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Flux function for multi-phase flow
 */ 
class bubbly_ff : public AdvDifFFWEnth {
private:
  int nelprops,nel,ndof,ndim,k_indx,e_indx,vg_indx,vl_indx,
    vl_indxe,vg_indxe;
  FastMat2 U,v_l,v_g,v_mix,Cp,Ajac,Id,Amoml,Amomg,Y,
    Djac,tmp1,Cjac,tmp2,tmp3,grad_v_l,strain_rate_l,
    grad_v_g,strain_rate_g,grad_k,grad_e,IdId,G_body,
    uintri,svec,tmp9,W_N,grad_alpha_l;
  FastMat2 Cpc,Ajacc,Djacc,Cjacc;
  double alpha_l,alpha_g,arho_l,arho_g,p,k,eps,
    visco_l,visco_g,visco_l_eff,visco_g_eff,
    C_mu,C_1,C_2,sigma_k,sigma_e,P_k,tau_fac,
    visco_t,temporal_stability_factor;

  double rho_l,rho_g;
    
public:
  bubbly_ff(NewElemset *elemset_);

  ~bubbly_ff();
  
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

};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// The elemset corresponding to the `bubbly_ff' flux function. 
class bubbly : public NewAdvDif {
public:
  /** Constructor, creates the fluc function object. 
      fixme:= should destroy the flux functin. 
  */
  bubbly() :  NewAdvDif(new bubbly_ff(this)) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class bubbly_bcconv : public NewBcconv {
public:
  bubbly_bcconv() : NewBcconv(new bubbly_ff(this)) {};
};

#endif

// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: aquifer.h,v 1.2 2002/01/15 16:25:54 mstorti Exp $
#ifndef AQUIFER_H
#define AQUIFER_H

#include "advective.h"

class NewAdvDifFFEnth : public NewAdvDifFF {
public:
  /** Allows updating the data for the object. 
      @param e (input) cofficients for updating the object
  */ 
  void update(const double *e) {};
  /** Computes the enthalpy vector from the state vector
      @param H (output) the enthalpy content vector
      @param U (input) the state vector
  */ 
  virtual void enthalpy(FastMat2 &H, FastMat2 &U)=0;
  /** Computes the product #(W_Cp_N)_(p,mu,q,nu) = W_p N_q Cp_(mu,nu)#
      @param W_Cp_N (output) size #nel# x #ndof# x #nel# x #ndof#
      @param W (input) weight function, size #nel#
      @param N (input) interpolation function, size #nel#
      @param w (input) scalar weight
  */ 
  virtual void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
			   double w)=0;
  /** Computes the product #(P_Cp)_(mu,nu) = (P_supg)_(mu,lambda) 
      Cp_(lambda,nu)#
      @param P_Cp (output) size #ndof# x #ndof#
      @param P_supg (input) matricial weight function, size #ndof# x #ndof#
  */ 
  virtual void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg)=0;
};

class IncludedEnthalpyFun : EnthalpyFun {
   NewAdvDifFFEnth 
}


class aquifer_ff : public NewAdvDifFFEnth {
  // K:= eta:=  
  /**  Basic properties, vertical position of aquifer bottom, hydraulic
       permeability, storativity 
  */
  Property eta_pr,K_pr,S_pr;
  /// Values of properties at an element
  double eta,K,S;
  /// Dimension of the problem (should be always 2)
  int ndim;
  /** Number of nodes per element, number of dof's per node, number of
      propoerties per element
  */
  int nel,ndof,nelprops;
public:
  /** This is called before any other in a loop and may help in
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
      fixme:= more doc here ...
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
  //@}

  /** @name Enthaply related functios */
  //@{
  /** Computes the enthalpy vector from the state vector
      @param H (output) the enthalpy content vector
      @param U (input) the state vector
  */ 
  void enthalpy(FastMat2 &H, FastMat2 &U);
  /** Computes the product #(W_Cp_N)_(p,mu,q,nu) = W_p N_q Cp_(mu,nu)#
      @param W_Cp_N (output) size #nel# x #ndof# x #nel# x #ndof#
      @param W (input) weight function, size #nel#
      @param N (input) interpolation function, size #nel#
      @param w (input) scalar weight
  */ 
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
			   double w);
  /** Computes the product #(P_Cp)_(mu,nu) = (P_supg)_(mu,lambda) 
      Cp_(lambda,nu)#
      @param P_Cp (output) size #ndof# x #ndof#
      @param P_supg (input) matricial weight function, size #ndof# x #ndof#
  */ 
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg);

  aquifer_ff(const NewAdvDif *e) : NewAdvDifFF(e) {}
};

class aquifer : public NewAdvDif {
 public:
  aquifer() :  NewAdvDif(new aquifer_ff(this)) {};
  // ~aquifer() { delete adv_diff_ff; }
};

#endif

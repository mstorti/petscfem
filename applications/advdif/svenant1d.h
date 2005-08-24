// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
// $Id: svenant1d.h,v 1.1 2005/08/24 01:52:49 mstorti Exp $
#ifndef PETSCFEM_SVENANT1D_H
#define PETSCFEM_SVENANT1D_H

#include "./advective.h"
#include "./advabso.h"
#include "./stream.h"
#include "./nonlres.h"

#include "nwadvdifj.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
/** Flux function for compressible viscous flow
 */
class svenant1d_ff : public AdvDifFFWEnth {
private:
  double Q,A,F,velmod,h_supg,
    w,wl,gravity,u,h;
  int ndim, nel, ndof, nelprops;
  const Elemset *old_elemset;
  FastMat2 U,W_N,Cp,Ajac;
public:
  svenant1d_ff(NewElemset *elemset=NULL);

  svenant1d_ff(Elemset *elemset=NULL);

  ~svenant1d_ff();

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

  void get_Cp(FastMat2 &Cp_a);

  void get_Ajac(FastMat2 &Ajac_a);

  void get_C(FastMat2 &C_a);

#define USE_COMP_P_SUPG
#ifdef USE_COMP_P_SUPG
  void comp_P_supg(FastMat2 &P_supg);
#endif

  /** The dimension of the corresponding `volume'
      elemenset minus ons. */
  int dim() const { return ndim; }

  /** Returns the Riemann Invariants and jacobians for Adv-Diff
      absorbent condition.  From #nel# nodes, nodes 1 to #nel-2# are
      ```real'', i.e. the nodes at the boundary, the #nel-1# node is
      reserved for the lagrange multipliers, and node #nel# is the
      nodes from where to take the reference state.
      @param U (input) (size #ndof#) state vector
      @param normal (input) (size #ndim#) normal to output surface.
      @param Rie (output) (size #ndof#) Riemman invariants for state #U#
      @param dRdU (output) (size #ndof x ndof#) Jacobian of Riemman invariants
      w.r.t. #U#
      @param (input) C (output) (size #ndof#) speeds for
      each Riemman characteristic variable. */
  void Riemann_Inv(const FastMat2 &U, const FastMat2 &normal,
		   FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C_);

  void
  compute_shock_cap_aniso(double &delta_aniso,
			  FastMat2 &jvec);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
/// The elemset corresponding to the `svenant1d_ff' flux function.
class svenant1d : public NewAdvDif {
public:
  /** Constructor, creates the flux function object.
      fixme:= should destroy the flux functin.
  */
  svenant1d() :  NewAdvDif(new svenant1d_ff(this)) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
class svenant1d_bcconv : public NewBcconv {
public:
  svenant1d_bcconv() : NewBcconv(new svenant1d_ff(this)) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
class svenant1d_abso : public AdvDiff_Abs_Nl_Res {
public:
  svenant1d_abso() :  AdvDiff_Abs_Nl_Res(new svenant1d_ff(this)) { }
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
class svenant1d_abso2 : public AdvectiveAbso {
public:
  svenant1d_abso2()
    :  AdvectiveAbso(new svenant1d_ff(this)) { }
};

#endif

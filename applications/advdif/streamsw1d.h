// -*- mode: C++ -*- 
// $Id: streamsw1d.h,v 1.3 2004/07/07 16:45:37 rodrigop Exp $
#ifndef STREAMSW1D_H
#define STREAMSW1D_H

#include "advective.h"
#include "stream.h"
#include "nonlres.h"
#include "godunov.h"

#define GETOPT_PROP(type,name,default) elemset->get_prop(name##_prop,#name)

//#define USE_COMP_P_SUPG

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The flux function for flow in a channel with arbitrary shape and
    using the shallow water 1d model
*/ 
class streamsw1d_ff : public AdvDifFFWEnth {
  // problem dimension
  int ndim;
  // element dim
  int ndimel;
  // dof=2 sw1d primitive variables (u h) 
  int ndof;
  // they must be here (channelshapes and frictionlaws)
  double area,wl_width,perimeter;
  // gravity
  double gravity;
  // allow scale stabilization term
  double  tau_fac;
  // scale friction term
  double cfric;
  // Threshold for h
  double h_min;
  // Threshold for v
  double vel_min;
  // Flux advec Jacobians
  FastMat2 A_jac;
  // Source Jacobian for newton loop
  FastMat2 Sf_jac;
  // Source friction terms Jacobian for newton loop
  FastMat2 C_jac;
  // Velocity vector in intrinsec coordinates
  FastMat2 Uintri;
  // Temp velocity vector
  FastMat2 UU;
  // Mass flux matrix
  FastMat2 flux_mass;
  // Momentum flux matrix
  FastMat2 flux_mom;
  // Cp enthalpy jacobian matrix
  FastMat2 Cp;
  // W_N  matrix (prod(W,N))
  FastMat2 W_N;
  // Temp matrix for flux functions
  FastMat2 tmp61,tmp2,tmp3;
  // Element iterator for hook
  ElementIterator elem_it;
  // Type of friction law used
  FrictionLaw *friction_law;
  // Pointer to the channel shape object
  ChannelShape *channel;

  //#define USE_A_JAC_DUMMY
#ifdef USE_A_JAC_DUMMY
  FastMat2 A_jac_dummy;
#endif

  double tmp_mask, adv_mask;

public:
  // Constructor
  streamsw1d_ff(const NewElemset *e);
  // Destructor
  ~streamsw1d_ff();
  /*
    This is called before any other in a loop and may help in
    optimization 
    @param options (input/output) this is used by the flux
    function writer for returning some options. Currently the only
    option used is #SCALAR_TAU#. This options tells the elemset
    whether the flux function returns a scalar or matrix
    #tau_supg#. 
  */ 
  void start_chunk(int &options);
  
  /* This is called before entering the Gauss points loop and may
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
  void comp_A_n(FastMat2 &A_jac_n, FastMat2 &normal);

  /** Computes fluxes, upwind parameters etc...
      fixme:= include more doc here ...
  */ 
  void compute_flux(COMPUTE_FLUX_ARGS);

  /** Computes the product #(N_P_C)_(mu,q,nu) 
      = w (P_supg)_(mu,lambda) C_(lambda,nu) N_q #
      @param N_P_C (output) size  #ndof# x #nel# x #ndof# 
      @param P_supg (input) SUPG perturbation function size #ndof# x #ndof#
      @param N (input) FEM interpolation function size #nel#
      @param w (input) a scalar coefficient
  */ 
  void comp_N_P_C(FastMat2 &N_P_C, FastMat2 &P_supg,
		  FastMat2 &N,double w);

  void comp_N_N_C(FastMat2 &N_N_C,FastMat2 &N,double w);

  //  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

  // Should be used only for Absorbing boundary conditions
  virtual void comp_A_jac_n(FastMat2 &A_jac_n, FastMat2 &normal) {
    assert(0);
  }

  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & grad_N,double w);

  void enthalpy(FastMat2 &H);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
  /* calcula el prod de la funcion de peso W con el jacobiano
     de la funcion de entalpia respecto de las variables primitivas y 
     la funcion de forma N en cada punto de gauss.
  */
  void comp_W_Cp_N(FastMat2 &W_Cp_N,const FastMat2 &W,const FastMat2 &N,
		   double w);
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:

  void comp_P_Cp(FastMat2 &P_Cp,const FastMat2 &P_supg);

  /** This stream elemset is essentially 1D.
      @return the dimension of this element that is 1
  */
  int dim() const { return 1; } //lo mismo que para el elemset stream de KWM

  /*
    Riemann Invariants calculus for Absorbent boundary conditions
    Rie, drdU  matrix declaration in boundary element routine rank=1 x ndof
  */
  void Riemann_Inv(const FastMat2 &U, const FastMat2 &normaln,
		   FastMat2 &Rie, FastMat2 &drdU, FastMat2 &C_U);

#ifdef USE_COMP_P_SUPG
  void comp_P_supg(FastMat2 &P_supg, FastMat2 &grad_N, FastMat2 &tau_supg) {
    P_supg.prod(A_jac,grad_N,-1,2,-1,1)
      .scale(tau_supg);
  }
#endif
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// The `streamsw1d' shallow water river or channel element.
class streamsw1d : public NewAdvDif {
public:
  streamsw1d() :  NewAdvDif(new streamsw1d_ff(this)) {
    // printf("En streamsw1d(): adv_diff_ff: %p\n",adv_diff_ff);
  } //constructor
};

class streamsw1d_abso : public AdvDiff_Abs_Nl_Res {
public:
  //streamsw1d_abso has n nodes [U_N, U_{N-1},U_{N-2}, .. , U_dummy]^t
  streamsw1d_abso() :  AdvDiff_Abs_Nl_Res(new streamsw1d_ff(this)/*,new streamsw1d(this)*/) {
    // printf("En streamsw1d_abso(): adv_diff_ff: %p\n",adv_diff_ff);
  } //constructor
};
class streamsw1d_abso_godunov : public AdvDiff_Godunov {
public:
  streamsw1d_abso_godunov() :  AdvDiff_Godunov(new streamsw1d_ff(this)) {
  } //constructor
};


#endif

// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: stream.h,v 1.3 2002/02/02 15:01:46 mstorti Exp $
#ifndef STREAM_H
#define STREAM_H

#include "advective.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AdvDifFFWEnth;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DummyEnthalpyFun : public EnthalpyFun {
  AdvDifFFWEnth *s;
public:
  DummyEnthalpyFun(AdvDifFFWEnth *s_) : s(s_) {}
  void enthalpy(FastMat2 &H, FastMat2 &U) ;
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class AdvDifFFWEnth : public NewAdvDifFF {
  DummyEnthalpyFun ef;
public:
  AdvDifFFWEnth(const NewElemset *elemset_=NULL) 
    : NewAdvDifFF(elemset_), ef(this) { enthalpy_fun = &ef; }
  virtual void enthalpy(FastMat2 &H, FastMat2 &U)=0;
  virtual void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w)=0;
  virtual void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg)=0;
};

/// Defines the shape of the channel
class ChannelShape {
  NewElemset *elemset;
  ChannelShape(NewElemset *e) : elemset(e) {}
public:
  ChannelShape *factory(NewElemset *e);
  
  /// Initialize properties perhap from the elemset table 
  virtual void init() {}

  /// Read local element properties
  virtual void 
  element_hook(NewElemset *elemset,ElementIterator element) {}
  
  /** For a given water depth (with respect to the bottom of the
      channel) give the fluid area, cross sectional water-line and wet
      channel perimeter.
      @param u (input) the water depth
      @param A (ouput) the fluid cross sectional area
      @param w (ouput) the fluid cross sectional water line
      @param P (ouput) the fluid cross sectional perimeter
  */
  virtual void area(double u,double &A,double &w,double &P)=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Rectangular shaped channel
class rect_channel : public ChannelShape {
  Property width_prop;
  double width;
public:
  void init();

  /// Read local element properties
  void element_hook(NewElemset *elemset,ElementIterator element) {}
  
  /** For a given water depth (with respect to the bottom of the
      channel) give the fluid area, cross sectional water-line and wet
      channel perimeter.
      @param u (input) the water depth
      @param A (ouput) the fluid cross sectional area
      @param w (ouput) the fluid cross sectional water line
      @param P (ouput) the fluid cross sectional perimeter
  */
  void geom_props(double u,double &A,double &w,double &P);
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class stream_ff : public AdvDifFFWEnth {
  double v;
  /// Pointer to the channel shape object
  ChannelShape *channel;
public:

  stream_ff(const NewAdvDif *e) 
    : AdvDifFFWEnth(e) }

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

  void enthalpy(FastMat2 &H, FastMat2 &U);
  void comp_W_Cp_N(FastMat2 &W_Cp_N,FastMat2 &W,FastMat2 &N,
		   double w);
  void comp_P_Cp(FastMat2 &P_Cp,FastMat2 &P_supg);

  /** This stream elemset is essentially 1D.
      @return the dimension of the advective elemset.
  */
  int dim() const { return 1; }
};

class rect_chan : public 

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class stream : public NewAdvDif {
 public:
  stream() :  NewAdvDif(new stream_ff(this)) {};
};

#endif

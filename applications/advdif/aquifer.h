// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: aquifer.h,v 1.9 2002/09/02 16:14:02 mstorti Exp $
#ifndef AQUIFER_H
#define AQUIFER_H

#include "advective.h"
#include "diff.h"

/** This is the flux function for a 
    quasi-harmonic equation with a conductivity proportional to the 
    difference between the free surface of the aquifer and its bottom 
    (a known quantity dependent on coordinates)
*/
class aquifer_ff : public DiffFF {
  // K:= eta:=  
  /**  Basic properties, vertical position of aquifer bottom, hydraulic
       permeability, storativity 
  */
  Property eta_pr,K_pr,S_pr;
  /// Values of properties at an element
  double phi,eta,K,S,rain,wet_aquifer_width_min;
  int dry_aquifer_stop;
  /// Dimension of the problem (should be always 2)
  int ndim;
  /** Number of nodes per element, number of dof's per node, number of
      propoerties per element
  */
  int nel,ndof,nelprops;
  /// Auxiliary matrices
  FastMat2 tmp,tmp1,tmp2,tmp3;

public:
  /** This is called before any other in a loop and may help in
      optimization 
  */ 
  void start_chunk();

  /** This is called before entering the Gauss points loop and may
      help in optimization. 
      @param element (input) an iterator on the elemlist. 
  */ 
  void element_hook(ElementIterator &element);

  /** This is called before each Gauss point.
      @param ipg (input) the Gauss point number (may be used to report errors)
  */ 
  void gp_hook(int ipg,const FastMat2 &U,const FastMat2 &grad_U);

  /** Computes fluxes, upwind parameters etc...
      fixme:= more doc here ...
  */ 
  void compute_flux(const FastMat2 &U,const FastMat2 &grad_U,
		    FastMat2 &fluxd,FastMat2 &G,
		    FastMat2 &H, FastMat2 &grad_H);

  /** Computes the product #(grad_N_D_grad_N)_(p,mu,q,nu) 
      = D_(i,j,mu,nu) (grad_N)_(i,p) (grad_N)_(j,q)#
      @param grad_N_D_grad_N (output) size #nel# x #ndof# x #nel# x #ndof# 
      @param grad_N (input) size #nel# x #ndof#
  */ 
  void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
			    FastMat2 & dshapex,double w);

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
  void comp_N_Cp_N(FastMat2 &N_Cp_N,FastMat2 &N, double w);

  /// This is called after all elements in a chunk
  virtual void end_chunk();

  aquifer_ff(const Diff *e) : DiffFF(e) {}
};

/// This is the elenset derived from the #aquifer_ff# flux function
class aquifer : public Diff {
 public:
  /// Constructor, creates the flux function object. 
  aquifer() :  Diff(new aquifer_ff(this)) {};
};

#endif

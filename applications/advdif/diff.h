// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: diff.h,v 1.5 2002/04/26 20:24:18 mstorti Exp $
#ifndef DIFF_H
#define DIFF_H

class Diff;

/// This is the flux function for a given physical problem. 
class DiffFF {
 public:
  /// The elemset associated with the flux function
  const Diff *elemset;
  /// Constructor from the elemset
  DiffFF(const Diff *elemset_=NULL) : elemset(elemset_) {};

  /// This flags whether a given task is processed or not
  int ask(const char *jobinfo,int &skip_elemset);

  /** This is called before any other in a loop and may help in
      optimization 
  */ 
  virtual void start_chunk() { }

  /** This is called before entering the Gauss points loop and may
      help in optimization. 
      @param element (input) an iterator on the elemlist. 
  */ 
  virtual void element_hook(ElementIterator &element) { }

  /** This is called before each Gauss point.
      @param ipg (input) the Gauss point number (may be used to report errors)
  */ 
  virtual void gp_hook(int ipg,const FastMat2 &U,const FastMat2 &grad_U) { }

  /** Computes fluxes, upwind parameters etc...
      fixme:= more doc here ...
  */ 
  virtual void compute_flux(const FastMat2 &U,const FastMat2 &grad_U,
			    FastMat2 &fluxd,FastMat2 &G,FastMat2 &H,
			    FastMat2 &grad_H)=0;

  /** Computes the product #(grad_N_D_grad_N)_(p,mu,q,nu) 
      = D_(i,j,mu,nu) (grad_N)_(i,p) (grad_N)_(j,q)#
      @param grad_N_D_grad_N (output) size #nel# x #ndof# x #nel# x #ndof# 
      @param grad_N (input) size #nel# x #ndof#
  */ 
  virtual void comp_grad_N_D_grad_N(FastMat2 &grad_N_D_grad_N,
				    FastMat2 & dshapex,double w) =0;

  /** Computes the enthalpy vector from the state vector
      @param H (output) the enthalpy content vector
      @param U (input) the state vector
  */ 
  virtual void enthalpy(FastMat2 &H, FastMat2 &U)=0;
  /** Computes the product #(W_Cp_N)_(p,mu,q,nu) = W_p N_q Cp_(mu,nu)#
      @param W_Cp_N (output) size #nel# x #ndof# x #nel# x #ndof#
      @param N (input) interpolation function, size #nel#
      @param w (input) scalar weight
  */ 
  virtual void comp_N_Cp_N(FastMat2 &N_Cp_N,FastMat2 &N,double w)=0;

  virtual ~DiffFF()=0;
};

/** Generic elemset for advective diffusive problems. 
    Several physical problems may be solved by defining the
    corresponding flux function object. 
*/
class Diff : public NewElemset { 
  /** A pointer to the flux function. Different
      physics are obtained by redefining the flux function
  */
  DiffFF *diff_ff;
public:
  /// Contructor from the pointer to the fux function
  Diff(DiffFF *diff_ff_=NULL) : diff_ff(diff_ff_) {};
  /// Destructor. Destroys the flux function object. 
  ~Diff() {delete diff_ff;}
  /// The assemble function for the elemset. 
  NewAssembleFunction new_assemble;
  /// The ask function for the elemset. 
  ASK_FUNCTION;
};



#endif

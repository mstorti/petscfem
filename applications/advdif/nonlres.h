// -*- mode: C++ -*-

#ifndef NONLRES_H  
#define NONLRES_H

#include "../../src/fem.h"
#include "advective.h"
#include <src/cloud.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:
class NonLinearRes : public NewElemset {
 public:
  ASK_FUNCTION;
  //  NewAssembleFunction new_assemble;
  /// Number of restrictions
  virtual int nres()=0;
  /// Initialize the elemset (maybe reads hash table)
  virtual void init()=0;
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. 
      @param jr (input) Number of restriction
      @param node (output) number of node for multiplier
      @param dof (output) number of field for multiplier   */ 
  virtual void lag_mul_dof(int jr,int &node,int &dof)=0;

  /** Computes the residual and jacobian of the function to be
      imposed. Usually you derive #NonLinearRes# and instantiate this
      function that defines the restriction to be imposed. 
      @param k (input) element number in elemset
      @param U (input) state vector at both nodes
      @param r (output) a vector of length #nres# containing the
      residuals
      @param lambda (input) the state of the Lagrange multipliers 
      @param jac (output) the jacobian of the residuals with respect
      of the node state. (size #nres x ndof#)
  */ 
  virtual void res(ElementIterator &element,FastMat2 &U,FastMat2 & r,
		   FastMat2 & lambda,FastMat2 & jac)=0;
  /// Make it pure virtual. 
  virtual ~NonLinearRes()=0;

};
 
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more that one node. 
*/ 
class LagrangeMult : public NewElemset {
 public:
  ASK_FUNCTION;
  //NewAssembleFunction new_assemble;
  /** Returns data (to be derived)
      @param nr (output) number of restrictions
      @param nfic (output) number of fictitious nodes
   */
  virtual int nres()=0;
  /** Return the node/dof pair to be used as lagrange multiplier for
      the #jr#-th restriction. 
      @param jr (input) Number of restriction
      @param node (output) number of node for multiplier
      @param dof (output) number of field for multiplier
  */ 
  virtual void lag_mul_dof(int jr,int &node,int &dof)=0;
  /// Initialize the elemset (maybe reads hash table)
  virtual void init()=0;
  /** Computes the residual and jacobian of the function to be
      imposed. Usually you derive #NonLinearRes# and instantiate this
      function that defines the restriction to be imposed. 
      @param element (input) element (iterator) in elemset
      @param U (input) state vector at all nodes
      @param r (output) a vector of length #nres*nel/2# containing the
      residuals for each restriction at each node.
      @param lambda (w) (input) the vector of reactions of the Lagrange multipliers 
      @param jac (output) the jacobian of the residuals with respect
      to the node state. (size #nel/2 * nres* nel/2 x ndof#)
  */ 
  virtual void res(ElementIterator &element,FastMat2 &U,FastMat2 & r,
		   FastMat2 & lambda,FastMat2 & jac)=0;
  /// Make it pure virtual. 
  virtual ~LagrangeMult()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This elemset provides the restriction (via Lagrange multipliers)
    of the nonlinear Dirichlet type b.c. at the wall of the form $k_w
    = u_*^2/C_\mu$ and similarly for $\epsilon_w$. 
*/
// classe para fijar caracteriticas de referencia para el prob de abs en adv-diff
class AdvDiff_Abs_Nl_Res : public NonLinearRes {
  //characteristics and Riemann invariant matrices
  FastMat2 RI_,C_U_,drdU_,drdU_ref,RI_ref,C_U_ref,
    xpe,cpe,RI_tmp,C_U_tmp,drdU_tmp,U_innodes,U_lagmul;
  /// Index (field numbers) of $u$ and $h$, number of dimensions
  int ndim,ndof,nu,nh,nel;
  /// Number of properties and restrictions
  int nprops,nr;
  double gravity;
  // u and h for Riemman Invariant of reference prop = u.dot(normal)+_(2(uh)^1/2)
  Property U_ref_prop ,normaln_prop;//normaln normal del nodo
  // u and h for Riemman Invariant of reference, inflow or outflow
  FastMat2 U_ref, normaln;
  // to get properties of elements
  ElementIterator element_m;
  // Poiter to adv-diff flux fun
  NewAdvDifFF *adv_diff_ff;
  // extrapolation cloud
  Cloud extr_cloud;
public:
  ASK_FUNCTION;
  NewAssembleFunction new_assemble;
  /// Number of restrictions (=2)
  int nres() {return ndof;}; 
  /// Initialize the elemset (maybe reads hash table)
  void init();
  // this returns the pair node/dof for lagrange multiplier for 
  // the #jr#-th restriction. 
  // @param jr (input) Number of restriction
  // @param node (output) number of node for multiplier
  // @param dof (output) number of field for multiplier
  void lag_mul_dof(int jr,int &node,int &dof);
  //element init
  void element_hook(ElementIterator &element);
  /// computes the residual and jacobian of the function to be imposed
  void res(ElementIterator &element,FastMat2 &U,FastMat2 &r,FastMat2 &lambda,
	   FastMat2 &jac);
  /// Constructor from the pointer to the flux function
  AdvDiff_Abs_Nl_Res(NewAdvDifFF *adv_diff_ff_=NULL) : 
    adv_diff_ff(adv_diff_ff_) {};
  /// Destructor
  ~AdvDiff_Abs_Nl_Res() {delete adv_diff_ff;};
};

#endif


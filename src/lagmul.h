// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: lagmul.h,v 1.6 2005/01/27 00:48:12 mstorti Exp $
#ifndef PETSCFEM_LAGMUL_H
#define PETSCFEM_LAGMUL_H

#define LagrangeMult GLagrangeMult

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element. 
    It may not work for restrictions that involve
    fields in more than one node. */ 
class LagrangeMult : public NewElemset {
 private:
  /// Stores internal matrix with coordinates of nodes
  FastMat2 xloc_m;
  /// Row dimension of coordinates vector
  int nu;
  /// The element actually visited
  int elem;
  /// Number of nodes per element
  int nel; 
  const Nodedata *nodedata_m;
 public:
  NewAssembleFunction new_assemble;
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
  /** Calls the #lm_initialize()# function for each elemset. */
  void initialize();
  /** Initialize the elemset. This is called in the
      LagrangeMult::initialize() function so that it is called before
      all chunks. And it is called even if there are not elements in
      this processor */
  virtual void lm_initialize() {}
  /** Initialize the elemset (maybe reads hash table). This is called before
      each element chunk. It is not called if there are not elements in
      this processor */
  virtual void init()=0;
  /** Computes the residual and jacobian of the function to be
      imposed. Usually you derive #NonLinearRes# and instantiate this
      function that defines the restriction to be imposed. 
      @param k (input) element number in elemset
      @param U (input) state vector at all nodes
      @param r (output) a vector of length #nres*nel/2# containing the
      residuals for each restriction at each node.
      @param w (input) the vector of reactions of the Lagrange multipliers 
      @param jac (output) the jacobian of the residuals with respect
      to the node state. (size #nel/2 * nres* nel/2 x ndof#)
  */ 
  virtual void res(int k,FastMat2 &U,FastMat2 & r,
		   FastMat2 & w,FastMat2 & jac)=0;
  /** Returns the coordinate of the nodes of the element. 
      @return a matrix with the coordinates of the nodes (size #nel*ndim#) */ 
  void get_xloc(ElementIterator element,
		FastMat2 &xloc,FastMat2 &H);
  /// Called after the loop over all elements
  virtual void close() {}
  /// Make it pure virtual. 
  virtual ~LagrangeMult()=0;

  virtual void 
  get_comp_flags(const char *jobinfo,
		 int &comp_mat,int &comp_mat_res)=0;
  virtual void
  get_data(arg_data_list &arg_data_v,
	   arg_data *&stateo,
	   arg_data *&staten,
	   arg_data *&retval,
	   arg_data *&retvalmat)=0;

  virtual void 
  get_data(arg_data_list &arg_data_v,
	   arg_data *&retvalmat)=0;

};

#undef LagrangeMult

#endif

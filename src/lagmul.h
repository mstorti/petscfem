// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: lagmul.h,v 1.12 2005/05/23 02:54:12 mstorti Exp $
#ifndef PETSCFEM_LAGMUL_H
#define PETSCFEM_LAGMUL_H

#include <src/penalize.h>

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
  ElementIterator element;
  /// Number of nodes per element
  int nel; 
  const Nodedata *nodedata_m;
  arg_data *stateo,*staten, 
    *retval, *retvalmat;
public:
  NewAssembleFunction new_assemble;
  /** Returns data (to be derived)
      @param nr (output) number of restrictions
      @param nfic (output) number of fictitious nodes
   */
  virtual int nres()=0;
  /** Return the node/dof pair to be used as lagrange
      multiplier for the #jr#-th restriction.
      @param jr (input) Number of restriction
      @param node (output) number of node for multiplier
      @param dof (output) number of field for multiplier
  */ 
  virtual void lag_mul_dof(int jr,int &node,int &dof)=0;
  /** Calls the #lm_initialize()# function for each
      elemset. */
  void initialize();
  /** Initialize the elemset. This is called in the
      LagrangeMult::initialize() function so that it is
      called before all chunks. And it is called even if
      there are not elements in this processor */
  virtual void lm_initialize() {}
  /** Initialize the elemset (maybe reads hash
      table). This is called before each element
      chunk. It is not called if there are not elements
      in this processor */
  virtual void init()=0;
  /** Computes the residual and jacobian of the function
      to be imposed. Usually you derive #NonLinearRes#
      and instantiate this function that defines the
      restriction to be imposed.
      @param k (input) element number in elemset
      @param U (input) state vector at all nodes
      @param r (output) a vector of length #nres*nel/2#
      containing the residuals for each restriction at
      each node.
      @param w (input) the vector of reactions of the
      Lagrange multipliers
      @param jac (output) the jacobian of the residuals
      with respect to the node state. (size #nel/2 *
      nres* nel/2 x ndof#)
  */ 
  virtual void res(int k,FastMat2 &U,FastMat2 & r,
		   FastMat2 & w,FastMat2 & jac)=0;
  virtual void set_ldf(FastMat2 &ldf_user,
		       vector<double> &ldf);
  /** Returns the coordinate of the nodes 
      of the element.       
      @param xloc (output) the coordinates of the nodes
      @param H (output) the vector of constant fields. */ 
  void get_xloc(FastMat2 &xloc,FastMat2 &H);
  /** Returns the old state at nodes. 
      @param Uold (output) the state of the nodes at the
      previous step. (size #nel*ndof#). */ 
  void get_old_state(FastMat2 &Uold);
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

  /// Initialize element
  virtual void element_hook(ElementIterator &element) {}
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** The `ROB' prefix means "Restriction-Object Based". 
    It means that the Lagrange multiplier object is based
    on a separate #Restriction# object, instead of deriving the
    the same Lagrange multiplier object. This allows the
    same restriction object to be used by several elemsets, 
    like #Penalize# and #LagrangeMult#. */ 
class ROBLagrangeMult : public LagrangeMult {
private:
  int nr;
  Restriction *restr;
  /// Number of nodes per element
  int nel, ndof, nelprops; 
  const Nodedata *nodedata_m;
  arg_data *stateo,*staten,*retval,*retvalmat;
public:
  int nres() { return nr; }
  void init() { 
    elem_params(nel,ndof,nelprops);
    nr = restr->init(nel,ndof,option_table(),name()); 
  }
  void lag_mul_dof(int jr,int &node,int &dof) {
    restr->lag_mul_dof(jr,node,dof);
  }
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac) {
    restr->res(k,U,r,w,jac);
  }
  void set_ldf(FastMat2 &ldf_user,
	       vector<double> &ldf) {
    restr->set_ldf(ldf_user,ldf);
  }
  void close() { restr->close(); }
  ROBLagrangeMult(Restriction *r=NULL) 
    : restr(r) { }
  ~ROBLagrangeMult() { 
    if (restr) {
      restr->close(); 
      delete restr; 
    }
  }
  // void element_hook(ElementIterator &element);
};

#undef LagrangeMult

#endif

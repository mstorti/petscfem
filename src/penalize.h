// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: penalize.h,v 1.1 2005/04/09 01:54:05 mstorti Exp $
#ifndef PETSCFEM_PENALIZE_H
#define PETSCFEM_PENALIZE_H

/** Enforces a restriction by adding a large term
    to the equation. The term is added proportional to
    a column-wise vector (as in the Lagrange multiplier
    formulation. */ 
class Restriction {
public:
  /** Initialize the object (maybe reads hash
      table). This is called before each element
      chunk. It is not called if there are not elements
      in this processor */
  virtual void init(TextHashTable *thash) { }
  /** Returns data (to be derived)
      @param nr (output) number of restrictions
      @param nfic (output) number of fictitious nodes
   */
  virtual int nres()=0;
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
  /// Called after the loop over all elements
  virtual void close() {}
  /// Make it pure virtual. 
  virtual ~Restriction()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Generic nonlinear restriction element, imposed
    via penalization. */ 
class Penalize : public NewElemset {
 private:
  /// Number of nodes per element
  int nel; 
  const Nodedata *nodedata_m;
  arg_data *stateo,*staten,*retval,*retvalmat;
  Restriction *restr;
 public:
  Penalize(Restriction *r=NULL) : restr(r) { }
  ~Penalize() { close(); if (restr) delete restr; }
  virtual void close()=0;
  NewAssembleFunction new_assemble;
  // virtual ~Penalize()=0;

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

#endif

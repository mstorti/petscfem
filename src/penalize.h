// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: penalize.h,v 1.8 2005/04/09 22:28:11 mstorti Exp $
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
      in this processor 
      @param nel (input) number of nodes per element
      @param ndof (input) number of fields per node
      @param thash (input) table of options
      @param name (input) the name of this object
      (probably equal to the element using it)
      @return number of restrictions to be applied
  */
  virtual int
  init(int nel,int ndof,
       TextHashTable *thash,const char *name) { }
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
/** Restriction loaded dynamically 
    (with #dlopen()# from a shared object file. */ 
class DLRestriction : public Restriction {
public:
  DLRestriction() { }
  ~DLRestriction() { }
  typedef 
  int InitFun(int nel,int ndof,
	      TextHashTable *thash,
	      const char *name,
	      void *&fun_data_a);
  typedef 
  void ResFun(int k,FastMat2 &U,FastMat2 & r,
	      FastMat2 & w,FastMat2 & jac,
	      void *fun_data_a);
  typedef void CloseFun(void *fun_data_a);
private:
  void *handle;
  void *fun_data;
  InitFun *init_fun;
  ResFun *res_fun;
  CloseFun *close_fun;
public:
  int init(int nel,int ndof,
	   TextHashTable *thash,const char *name);
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac) {
    (*res_fun)(k,U,r,w,jac,fun_data);
  }
  /// Called after the loop over all elements
  void close() { (*close_fun)(fun_data); }
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
  /// The element actually visited
  int elem;
  ElementIterator element;
 public:
  Penalize(Restriction *r=NULL) : restr(r) { }
  ~Penalize() { 
    if (restr) {
      restr->close(); 
      delete restr; 
    }
  }
  virtual void close() { }
  NewAssembleFunction new_assemble;

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

#define DL_GENERIC_RESTRICTION(prefix)			\
extern "C"						\
void prefix##_init_fun(int nel,int ndof,		\
		       TextHashTable *thash,		\
		       const char *name,		\
		       void *&fun_data) {		\
  fun_data = new prefix;				\
  ((prefix *)fun_data)->init(nel,ndof,thash,name);	\
}							\
							\
extern "C" void						\
prefix##_res_fun(int k,FastMat2 &U,FastMat2 & r,	\
	      FastMat2 & w,FastMat2 & jac,		\
	      void *fun_data) {				\
  ((prefix *)fun_data)->res(k,U,r,w,jac);		\
}							\
							\
extern "C" void						\
prefix##_close_fun( void *fun_data) {			\
  ((prefix *)fun_data)->close();			\
}

#endif

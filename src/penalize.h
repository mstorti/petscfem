// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: penalize.h,v 1.15 2005/08/05 01:54:50 mstorti Exp $
#ifndef PETSCFEM_PENALIZE_H
#define PETSCFEM_PENALIZE_H

/** Enforces a restriction by adding a large term
    to the equation. The term is added proportional to
    a column-wise vector (as in the Lagrange multiplier
    formulation). */ 
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
  /** Return the node/dof pair to be used as lagrange
      multiplier for the #jr#-th restriction.
      @param jr (input) Number of restriction
      @param node (output) number of node for multiplier
      @param dof (output) number of field for multiplier
  */ 
  virtual void lag_mul_dof(int jr,int &node,int &dof) { 
    assert(0); // force user to instantiate it
  }
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
  virtual void uold(const FastMat2 &Uold) { }
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
  void LagMulDofFun(int jr,int &node,
		    int &dof,void *fun_data_a);
  typedef
  void UoldFun(const FastMat2 &uold,void *fun_data_a);
  typedef 
  void ResFun(int k,FastMat2 &U,FastMat2 & r,
	      FastMat2 & w,FastMat2 & jac,
	      void *fun_data_a);
  typedef void CloseFun(void *fun_data_a);
private:
  void *handle;
  void *fun_data;
  InitFun *init_fun;
  LagMulDofFun *lag_mul_dof_fun;
  UoldFun *uold_fun;
  ResFun *res_fun;
  CloseFun *close_fun;
public:
  int init(int nel,int ndof,
	   TextHashTable *thash,const char *name);
  void lag_mul_dof(int jr,int &node,int &dof) {
    assert(lag_mul_dof_fun);
    (*lag_mul_dof_fun)(jr,node,dof,fun_data);
  }
  void uold(const FastMat2 &Uold) {
    (*uold_fun)(Uold,fun_data);
  }
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class DLBaseRestriction {
public:
  virtual 
  int init(int nel,int ndof_a,
	   TextHashTable *thash,const char *name) { return 0; }
  virtual
  void res(int k,FastMat2 &U,FastMat2 & r,
	   FastMat2 & w,FastMat2 & jac)=0;
  virtual
  void lag_mul_dof(int jr,int &node,int &dof) { assert(0); }
  virtual 
  void uold(const FastMat2 &Uold) { }
  virtual
  void close() { }
};

// Takes the generic pointer `fun_data' and converts it to
// the object.  WARNING: Can't use RTTI (dynamic_cast<>)
// here because `DLBaseRestriction' is not a derived class
// of `void'.  Using RTTI would ensure that the proper
// object is loaded, and also would let defining default
// behaviour for the objects.  `fun_data' should be declared
// as a base class for `DLBaseRestriction'.
#define dl_penal_convert_to_class		\
DLBaseRestriction *obj =			\
(DLBaseRestriction *)(fun_data);		\
assert(obj);

// Wraps the object interface to a C interface that
// can be dynamically laoded. Each C function calls the
// corresponding function on the object hided in the
// `fun_data' argument. 
#define DL_GENERIC_RESTRICTION(prefix)			\
extern "C"						\
void prefix##_init_fun(int nel,int ndof,		\
		       TextHashTable *thash,		\
		       const char *name,		\
		       void *&fun_data) {		\
  fun_data = new prefix;				\
  dl_penal_convert_to_class;				\
  obj->init(nel,ndof,thash,name);			\
}							\
							\
extern "C" void						\
prefix##_lag_mul_dof_fun(int jr,int &node,		\
			 int &dof,void *fun_data) {	\
  dl_penal_convert_to_class;				\
  obj->lag_mul_dof(jr,node,dof);			\
}							\
							\
extern "C" void						\
prefix##_uold_fun(const FastMat2 &uold,			\
		  void *fun_data) {			\
  dl_penal_convert_to_class;				\
  obj->uold(uold);					\
}							\
							\
extern "C" void						\
prefix##_res_fun(int k,FastMat2 &U,FastMat2 & r,	\
	      FastMat2 & w,FastMat2 & jac,		\
	      void *fun_data) {				\
  dl_penal_convert_to_class;				\
  obj->res(k,U,r,w,jac);				\
}							\
							\
extern "C" void						\
prefix##_close_fun( void *fun_data) {			\
  dl_penal_convert_to_class;				\
  obj->close();						\
}

#endif

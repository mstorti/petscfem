// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: elemset.h,v 1.21 2002/08/31 20:14:57 mstorti Exp $

#ifndef ELEMSET_H
#define ELEMSET_H

#include "libretto.h"
#include <glib.h>

#include "arglist.h"
#include "getprop.h"

enum ElemsetIteratorMode {
  ALL                  = 0x00001,
  INCLUDE_GHOST        = 0x00010,
  DO_NOT_INCLUDE_GHOST = 0x00100
};

class ElementList;
class ElementIterator;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Stores sets of elements with similar types and properties. 
    @author M. Storti
*/ 
class Elemset {
public:
  // type of element
  char *type;
  /// table of connectivities
  int *icone;
  ///number of elements in the elemset
  int nelem;
  ///number of nodes per element
  int nel;
  ///number of degrees of freedom per node
  int ndof;
  /// mesh partitioning as computed by Metis 
  int *epart;
  /// flag indicating whether this  is the fat elemset or not
  int isfat;
  /// number of elements in this processor
  int nelem_here;
  /// number of double properties in the per-element properties table
  int nelprops;
  /// number of integer properties in the per-element properties table
  int neliprops;
  /// number of additional double properties 
  int nelprops_add;
  /// number of additional integer properties
  int neliprops_add;
  /// double per-element properties table
  double *elemprops;
  /// int per-element properties table
  int *elemiprops;
  /// additional double per-element properties table
  double *elemprops_add;
  /// additional int per-element properties table
  int *elemiprops_add;
  /// This is a ``loccker'' for each element
  void **local_store;
  /// properties hash table 
  TextHashTable *thash;
  /// hash of properties in the per element double properties table
  GHashTable *elem_prop_names;
  /// hash of properties in the per element int properties table
  GHashTable *elem_iprop_names;

  /// Makes some initialization from the hash table
  virtual void initialize() {}

  /// Returns the list of "real" nodes (this is used for computing the graph)
  virtual int real_nodes(int iele,const int *&nodes);

  /** Assembles residuals, matrices, scalars and other quantities. 
      @author M. Storti
      @param retval (output) here returns contributions vectors
      (residuals), and matrices
      @param nodedata (input) vector with properties per node
      @param locst (input) state vectors localized to elements 
      @param locst2 (input) same for alternative  state vector (this may be
      used in temporal integration), etc...
      @param dofmap (input) maps the node/field representation to the vector
      state
      @param ijob (input) tells the global assemble function  what kind of job
      should be done (assemble vector, profile or matrix). To be
      processed by the global assemble function. 
      @param jobinfo (input) tells the routine element what kind of
      matrix or vector has to be assembled
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
  */
  virtual int assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
		       const char *jobinfo,int myrank,
		       int el_start,int el_last,int iter_mode,
		       const TimeData *time_data) {
    
    printf("assemble: not known Elemset\n"); exit(1);};

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Ask the elemset if it should be processed for this jobinfo. 
      @author M. Storti
      @param jobinfo (input) The name of the task.
      @param answer (output) 1 = process (0 = do not process) this
      elemset. 
  */ 
  virtual int ask(const char *jobinfo,int &skip_elemset) {
    
    // By default process the elemset. 
    skip_elemset = 0;
    return 0;
  };
    

  /// dynamic aray with ghost elements
  Darray *ghost_elems;
  /// print info of this elemset
  void print();
  
  /** Returns the amount of work needed to process this element.
      @return the weight of the processor
  */ 
  virtual double weight();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Builds localized vectors (in node/field) representation of state
      vectors. 
      
      @author M. Storti
      @param nel (input) number of nodes per element
      @param ndof (input) number of degrees of freedom per node
      @param dofmap (input) the dofmap of the mesh
      @param locst (output) localized vector to be assembled
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
      @param time_data (input, def=NULL) an external parameter in order to compute
      external boundary conditions, etc...
  */ 
  int download_vector(int nel,int ndof,Dofmap *dofmap,
		      arg_data &argd,
		      int myrank,int el_start,int el_last,int iter_mode,
		      const TimeData *time_data=NULL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Put localized values (nod/field representation) returned by
      elemset assembles in the global residual vector. 
      
      @author M. Storti
      @param nel (input) number of nodes per element
      @param ndof (input) number of degrees of freedom per node
      @param dofmap (input) the dofmap of the mesh
      @param retval (input) localized values (nod/field
      representation) returned by elemset assembles 
      @param A (input/output) assemble matrix contributions on it (if
      ijob/jobinfo indicates so). 
      @param da (output) assemble matrix profile on this (if
      ijob/jobinfo indicates so). 
      @param vec (output) assemble vector contributions on it (if
      ijob/jobinfo indicates so).
      @param myrank (input) identifies this processor
      @param el_start (input) low end of the range of elements to be
      processed
      @param el_last (input) high end of the range of elements to be
      processed
      @param iter_mode (input) include or not ghost elements
      @param klocc (input) if computing difference finite jacobian, is
      the local node perturbed 
      @param klocc (input) if computing difference finite jacobian, is
      the local d.o.f. perturbed 
  */
  int upload_vector(int nel,int ndof,Dofmap *dofmap,
		    int options,arg_data &argd,int myrank,
		    int el_start,int el_last,int iter_mode,
		    int klocc=0,int kdofc=0);

  /** Return the ``locker'' for the element, that is a (void *) where to
      put data. Only for elements local to this processor. 
      @author M. Storti
      @param global_elem (input) the index of the element in the elemset
      @return the address (void *) of the locker
  */
  void *& local_store_address(int global_elem) {
    return local_store[global_elem];
  }

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns coordinates and node data for the element
      @author M. Storti
      @param element (input) iterator to the element for which 
      to return the data
      @param nodedata (input) the node coordinate info
      @param xloc (input) the coordinates of the element nodes
      @param Hloc (input) the auxiliary data for the nodes
  */ 
  void element_node_data(const ElementIterator &element,
			 const Nodedata *nodedata,
			 double *xloc,double *Hloc) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns element properties
      @author M. Storti
      @param element (input) iterator to the element for which 
      to return the data
      @return pointer to an array of nelprops doubles for the element
  */ 
  double *
  element_props(const ElementIterator &element) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns node indices for the element
      @author M. Storti
      @param element (input) iterator to the element for which 
      to return the data
      @param connect (input) indices of the element nodes
  */ 
  void element_connect(const ElementIterator &element,
		       int *connect) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized vector values for a given element
      @param element (input) iterator to the element for which 
      to return the data
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  const double *
  element_vector_values(const ElementIterator &element,
			arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized element vector contributions for a given element
      @param element (input) iterator to the element for which 
      to return the data
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  element_ret_vector_values(const ElementIterator &element,
			    arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized element contributions to FD 
      jacobian for a given element
      @param element (input) iterator to the element for which 
      to return the data
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  element_ret_fdj_values(const ElementIterator &element,
			 arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized element contributions to FD 
      jacobian for a given element
      @param element (input) iterator to the element for which 
      to return the data
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  element_ret_mat_values(const ElementIterator &element,
			 arg_data &ad) const;

  /** Returns element values
      @author M. Storti
      @param nel_ (output) the number of nodes connected to an element
      @param ndof_ (output) the number of dofs for each node
  */
  void elem_params(int &nel_,int &ndof_,int &nelprops_) const {
    nel_=nel; ndof_=ndof; nelprops_ = nelprops;
  }

  const char *name();

  friend class ElementList;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  // The following will be declared private in the future
  //private: 

  /// Stores temporarily element connectivities
  int *elem_conne;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Global assemble function. Loops over all elemsets and assembles
    element contributions to amtrices, vectors and profiles. 
    Same as assemble, but for computine profiles.
    @author M. Storti
    @param mesh (input) the mesh to be processed
    @param argl (input/output) list of arguments, (vector matrices) on input
    and output. 
    @param dofmap (input) the dofmap of the mesh
    @param jobinfo (input) tells the routine element what kind of
    matrix or vector has to be assembled
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @return error code
*/
int assemble(Mesh *mesh,arg_list argl,Dofmap *dofmap,const char *jobinfo,
	     const TimeData *time_data=NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Const iterator for looping over elements in an elemset.
class ElementList {
private:
  /// Pointer to the elemset
  const Elemset *elemset;
  /// Start range of values while iterating
  int first;
  /// Ends range of values while iterating
  int last;
  /// Iteration mode
  ElemsetIteratorMode mode;
public:
  /// Constructor
  ElementList(const Elemset *elemset_,
	      int first_,int last_,ElemsetIteratorMode mode_) :
    elemset(elemset_), first(first_), last(last_), mode(mode_) {};
  /// Allow iterator to access this class
  friend class ElementIterator;
  /// Begin of chunk
  ElementIterator begin(void) const;
  /// End of chunk
  ElementIterator end(void) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Iterators for looping over elements in an elemset.
class ElementIterator {
private:
  /// The list to which it belongs
  const ElementList *elemlist;
  /// Rank of element in the chunk
  int rank_in_chunk;
  /// Rank of element in the elemset
  int rank_in_elemset;
public:
  /// Default Constructor
  ElementIterator() {};
  /// Constructor
  ElementIterator(const ElementList *el, const int rie,const int ric) 
    : elemlist(el), rank_in_elemset(rie), rank_in_chunk(ric) {};
  /// Prefix increment operator
  ElementIterator & operator++(void);
  /// Postfix Increment operator
  ElementIterator operator++(int);
  /// Not equal operator
  int operator!=(const ElementIterator &other) const;
  /// Return position in chunk
  void position(int &pos_in_elemset, int &pos_in_chunk) const {
    pos_in_elemset = rank_in_elemset; pos_in_chunk = rank_in_chunk;
  }
  /// Returns true/false whether the current element is
  /// in the restricted list or not
  int is_valid() const;
  /// Advances iterator to 
  /// in the restricted list or not
  void advance_to_valid();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns coordinates and node data for the element
      @author M. Storti
      @param nodedata (input) the node coordinate info
      @param xloc (input) the coordinates of the element nodes
      @param Hloc (input) the auxiliary data for the nodes
  */ 
  void node_data(const Nodedata *nodedata,double *xloc,double *Hloc);
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized vector values for a given element
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  const double *
  vector_values(arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized vector element contributions for a given element
      @param element (input) iterator to the element for which 
      to return the data
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  ret_vector_values(arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized vector element contributions to FD 
      jacobian for a given element
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  ret_fdj_values(arg_data &ad) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns localized matrix element contributions
      for a given element
      @param ad (input) the argument in the arglist.
      @return a pointer to the element values
      @author M. Storti
  */ 
  double *
  ret_mat_values(arg_data &ad) const;

  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns element properties
      @author M. Storti
      @return pointer to an array of nelprops doubles for the element
  */ 
  double * props();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** @name compute\_prof package */
//@{
/** nodes to store the profile of a matrix. 
    @author M. Storti
*/ 
class Node {
public:
  /// points to next dof connected to this one
  int next,
    /// the connected dof
    val;
  /// constructor
  Node(int next_=-1 ,int val_=0 ) : next(next_), val(val_) {};
  /// printer
  void print(void) {printf(" (%d,%d)\n",next,val);};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Inserts a node in the profile. 
    @author M. Storti
    @param da (input/output) dynamic array containing the profile
    @param j (input) row index
    @param k (output) column index
*/ 
void node_insert(Darray *da,int j,int k);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Compare function for sort and then binary search on the profile. 
    @author M. Storti
    @param left (input) pointer to left element to be compared
    @param right (input) pointer to right element to be compared
    @param args (not used) as required by libretto routines
*/ 
int int_cmp (const void *left,const void *right, void *args);
//@}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Flags if the current element is to be computed his contribution. 
    Basically it process all elements in the current
    processor. However, depending on iter\_mode it may include also
    ghot-elements. 

    @author M. Storti
    @param iele (input) element to be tested
    @param elemset (input) this elemset 
    @param myrank (input) this processor index
    @param iter\_mode (input) flags whether or not include
    ghost-elements. 
    @return boolean flag indicating whether process or not this element
*/ 
int compute_this_elem(const int & iele,const Elemset *elemset,const int & myrank,
		      int iter_mode);

#define ASSEMBLE_FUNCTION \
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata, \
	       Dofmap *dofmap,const char *jobinfo,int myrank, \
	       int el_start,int el_last,int iter_mode, \
	       const TimeData *)

#define ELEMSET_CLASS(name)			\
class name : public Elemset {			\
public: ASSEMBLE_FUNCTION;			\
}

#define ASK_FUNCTION int ask(const char *jobinfo,int &answer)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Loops many times executing assmble() with the same arguments and
    issuing a performance indication. 
    @author M. Storti
    @param mesh (input) the mesh to be processed
    @param argl (input/output) list of arguments, (vector matrices) on input
    and output. 
    @param dofmap (input) the dofmap of the mesh
    @param jobinfo (input) tells the routine element what kind of
    matrix or vector has to be assembled
    @param time_data (input, def=NULL) an external parameter in order to compute
    external boundary conditions, etc...
    @return error code
*/ 
int measure_performance_fun(Mesh *mesh,arg_list argl,
			    Dofmap *dofmap,const char *jobinfo,const TimeData
			    *time_data=NULL);

typedef void 
NewAssembleFunction(arg_data_list &arg_datav,const Nodedata *nodedata,
		    const Dofmap *dofmap,
		    const char *jobinfo,const ElementList &elemlist,
		    const TimeData *time_data);

/// The generic ElemProperty class
class Property {
public:
  int indx;
  vector<double> val;
  double *ptr;
  int length;
  Property() : indx(-1), length(0), ptr(NULL) {val.clear();};
};

#if 0
class Property {
public:
  int indx;
  double *val;
  Property() : indx(-1), val(NULL) {};
  ~Property() {if (val) delete[] val;};
};
#endif

/** This is an adaptor to the old Elemset class
    @author M. Storti
*/
class NewElemset : private Elemset {
  /// This is the adaptor to the old assemble function.
  int assemble(arg_data_list &arg_datav,Nodedata *nodedata,Dofmap *dofmap,
	       const char *jobinfo,int myrank,
	       int el_start,int el_last,int iter_mode,
	       const TimeData *time_data);

  /** @name The new get-prop methods */
  //@{
  /// The vector containing the double values
  vector<double> propel;
  /// The first position in propel
  double *begin_propel;
  /// The indices 
  vector<int> elprpsindx;

  //@}
public:
  /// The new assemble function
  virtual void 
  new_assemble(arg_data_list &arg_datav,const Nodedata *nodedata,const Dofmap *dofmap,
	       const char *jobinfo,const ElementList &elemlist,
	       const TimeData *time_data) {
    printf("assemble: not known New Elemset\n"); exit(1);
  };

  Elemset::initialize;
  Elemset::elem_params;
  Elemset::local_store_address;
  Elemset::element_node_data;
  Elemset::element_props;
  Elemset::element_connect;
  Elemset::element_vector_values;
  Elemset::element_ret_vector_values;
  Elemset::element_ret_fdj_values;
  Elemset::element_ret_mat_values;
  // Elemset::thash;
  void get_entry(const char *k,const char *&v) const {
    thash->get_entry(k,v);};

  int get_int(const char *name,
	      int &retval,int defval=0,int n=1) const {
    return ::get_int(thash,name,&retval,defval,n);
  };
  
  int get_double(const char *name,
		 double &retval,int defval=0,int n=1) const {
    return ::get_double(thash,name,&retval,defval,n);
  };

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a vector of double properties whose length is unknown. 
      @author M. Storti
      @param name (input) the name of the property
      @param retval (input) the returned vector
      @param defval (input) Controls the action if no entry is found
       in the hash table. If #defval==0#, let retval unchanged, so that
      you have to set it to its default value. Give an error if
      #defval!=0# 
   */ 
  int get_vec_double(const char *name,
		     vector<double> &retval,int defval=0) const;

  int get_string(const char *name,
		 string &ret,int defval=0,int n=1) const {
    return ::get_string(thash,name,ret,defval,n);
  };

  /// Creates a Property object from his name
  void get_prop(Property &prop,const char *prop_name,int n=1) const;
  double prop_val(ElementIterator &element,Property &prop) const;
  const double *prop_array(ElementIterator &element,Property &prop) const;
};

#if 0
#define NEW_ASSEMBLE_FUNCTION \
  int assemble(arg_data_list &arg_data_v,Nodedata *nodedata, \
	       Dofmap *dofmap,const char *jobinfo,int myrank, \
	       int el_start,int el_last,int iter_mode, \
	       const TimeData *)
#endif

#endif

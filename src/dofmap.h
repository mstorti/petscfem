// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: dofmap.h,v 1.24.10.1 2007/02/19 20:23:56 mstorti Exp $
 
#ifndef DOFMAP_H
#define DOFMAP_H

#include <vector>
#include <map>
#include <string>
#include <iostream>

#include <petsc.h>
#include <petscvec.h>

#include <src/part.h>
#include <src/texthash.h>
#include <src/fstack.h>
#include <src/idmap.h>
#include <src/fastlib.h>
#include <src/fastlib2.h>
#include <src/dvector.h>

class TimeData {};

// class TimeData {
// public:
//   virtual double time() const;
// };

class Time : public TimeData{public:
  Time(double t=0.0) : time_(t) {};
  Time(const Time& t) : time_(t.time_) {};
  double time() const {return time_;};
  void set(const double t)  {time_=t;};
  void inc(const double dt) {time_+=dt;};
  operator double() const {return time_;};
  Time& operator=(const Time& t) { time_=t.time_; return *this;}
private:
  double time_;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** An amplitude is basically a function object that returns the
    Dirichlet value for a given time. For instance function `sine'  
    with constant parameters  `omega' and `phase', depending on time. 
    @author M. Storti
*/ 
class Amplitude {
public:
  virtual ~Amplitude() {}
  static Amplitude *factory(char *& label,TextHashTable *tht_=NULL);
  static Amplitude *old_factory(char *& label,FileStack &fstack);
  /// Eval the amplitude of the function at this time. 
  virtual double eval(const TimeData *time_data);
  /// Eval the amplitude of the function at this time (needs node and field)
  virtual double eval(const TimeData *time_data,int node,int field);
  /** Callback function defined by the user -- returns
      whether this amplitude function needs to be passed the
      node/field combination. 
  */
  virtual int needs_dof_field_q();
  /** Initializes the object. Table t should be deleted if not
      incorporated in the created object. If it is deleted set it to
      NULL. 
  */
  virtual void init(TextHashTable *t) { delete t; t=NULL; };
  /// prints the amplitude entry. 
  virtual void print(void) const=0; 
};

/** @name dofmap.
    This is not documented since it will suffer from major revision when
    using the idmap class
*/

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** This is tha basic element in the lists that represents a
    constraint.  The constraint is defined by passing to the dofmap a
    list of node/field coefficients, then the constraint is this
    linear combination set to 0. If a non-homogeneous constraint is
    needed then a fixed dof should be introduced. 
    @author M. Storti
*/ 
struct constraint_entry {
  /// Node number. 
  int node;
  /// Field number. 
  int field;
  /// Coefficient for this node/field. 
  double coef;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Constraints are set passing a list of node, field,
    coefficients. This class implements such lists. 
*/
class Constraint : public vector<constraint_entry>  {

public:
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Add an entry to the list. 
      @author M. Storti
      @param node (input) node number
      @param field (input) field number
      @param coef (input) coefficient in the linear combination 
   */ 
  void add_entry(int node,int field,double coef);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Empties the list of the linear combination. 
      @author M. Storti
   */ 
  void empty(); 
};

// forward declaration
class Dofmap;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Defines a fixation (typically Dirichlet b.c.'s).  The double is
    the "spatial" amplitude and the int is a pointer (or something)
    that describes the temporal evolution.
*/
//typedef pair<double,int> fixation_entry;
class fixation_entry {
public:
  friend class Dofmap;
private:
  double val;
  int edof;
  Amplitude *amp;
public:
  fixation_entry(double val_=0.,Amplitude *amp_=NULL,int edof_a=-1) :
    val(val_), edof(edof_a), amp(amp_) {};
  // double value(const TimeData *time_data) const;
  void print(void) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class q_entry {
public:
  int node,kdof,keq;
  double coef;
  q_entry(const int node_=0,const int kdof_=0,const int keq_=0,
	  const double coef_=0.)
    : node(node_), kdof(kdof_), keq(keq_), coef(coef_) {};
  void print();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int q_entry_cmp (const void *left,const void *right, void *args);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class fixa_entry {
public:
  int node,kdof;
  double val;
  fixa_entry(const int node_=0,const int kdof_=0,const double val_=0.)
    : node(node_), kdof(kdof_),val(val_) {};
  void print();
};

int fixa_entry_cmp (const void *left,const void *right, void *args);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class sp_entry {
public:
  int node,kdof,gdof;
  double coef;
  sp_entry (int node_=0,int kdof_=0, int gdof_=0,
	    double coef_=0.) : node(node_),kdof(kdof_),gdof(gdof_),coef(coef_) {};
};

//-------<*>-------<*>-------<*>-------<*>-------<*>------- 
/** Stores the mapping between the node/field representation and
    unknowns in MPI vectors and matrices.

    @author M. Storti
*/
class Dofmap : public DofPartitioner {

private:
  /** This maps a node/dof to an eq. or either an entry in the 
   non-regular list. */
  dvector<int> idmap2_dv;
  int *idmap2;
  dvector<int> special_ptr_dv;
  int *special_ptr;
  dvector<int> sp_eq_dv;
  int *sp_eq;
  dvector<double> coefs_dv;
  double *coefs;
  double one_coef;
  Vec x,z;
  vector<double> w;

public:
  /// MPI communicator (default: PETSCFEM_COMM_WORLD)
  MPI_Comm comm;

  /// number of nodes
  int nnod;

  /// number of equations
  int neq;

  /** start range of dof's for this processor. Dof's on this processor
      are dof1<= dof <=dof2. (dof in base 1). 
   */
  int dof1;

  /// end range of dof's for this processor. (See arg dof1). 
  int dof2;

  /// number of independent fixations (entries in the `fixed' table)
  int  neqf;

  /** total number of independent variables (fixed + free) = column
      dimension of Q
  */
  int neqtot;

  /// number of degrees of freedom per node
  int ndof;

  /** identifier of whether the node/field pair is free,
	fixed, or special. (This will be included in an
	idmap object in the future.) */
  // fixme:= borrar este despues de pasar a idmap!!
  int *ident;

  /// The STL vector containing ghost\_dofs (ordered, base 0)
  vector<int> *ghost_dofs;

  /// this will replace ident in a future
  idmap *id;
  /// points to fixations
  Darray *fixa;
  /// points to special rows
  Darray *q;
  /** stores fixations. `fixa' will change to this. 
      In the Q matrix, column indices greater than `neq' point to
      fixation j-neq-1
   */
  vector<fixation_entry> fixed;
  /// flags whether the dofmap is synchronized or not
  //  int synchro; // indicates if `q' and `fixa' are sorted
  /// startproc[myrank] id the first row index in processor `myrank'
  int *startproc;
  /// number of lines (unknowns) in processor `myrank'
  int *neqproc;
  /// number of processors
  int size;
  /// weight per processor
  float *tpwgts;
  /// node partitioning (as returned by Metis)
  int *npart;
  /// Scatter to convert to sequential vector with ghost values
  VecScatter ghost_scatter;
  /// scatter to print
  VecScatter scatter_print;

  // This is used temporarily in order to store the mappings between
  // edofs and fixations
  map<int,int> fixed_dofs;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a row of the mapping matrix, only free (jeq<=neq) dofs are returned.
      @author M. Storti
      @param node (input) the node number to get the row
      @param kdof (input) the field number to get the row
      @param row (output) the rwtrieved row */ 
  void get_row_free(int const & node,int const & kdof,row_t &row) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a row of the maping matrix.
      @author M. Storti
      @param node (input) the node number to get the row
      @param kdof (input) the field number to get the row
      @param row (output) the retrieved row */ 
  void get_row(int const & node,int const & kdof,row_t &row) const;

  void get_row(int const & node,int const & kdof,IdMapRow &row) const;

  /** Fast version. Uses an internal array. ADD ARG DOC ... */ 
  void get_row(int node,int kdof,int &ndof,const int **dof,
	       const double**coef) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets a row of the maping matrix.
      @author M. Storti
      @param node (input) the node number to get the row
      @param kdof (input) the field number to get the row
      @param row (input) the retrieved row */ 
  void row_set(const int & node,const int & kdof,const row_t &row);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a dof value from a state vector scattered to this processor.. 
      @author M. Storti
      @param jeq (input) the column index to be retrieved. May be free
      (1 < jeq <= neq) or fixed (neq < jeq <= neqtot)
      @param sstate (input) the state vector
      @param ghost_vals (input) double array containing scattered ghost values
      @param time_data (input) a pointer to a struct that tipically is a time
      @return the double value corresonding to that dof */ 
  double get_dofval(const int & jeq,double const *sstate, double const
		    *ghost_vals,const TimeData *time_data) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a dof value from state vector and scattered ghost values. 
      @author M. Storti
      @param jeq (input) the column index to be retrieved. May be free
      (1 < jeq <= neq) or fixed (neq < jeq <= neqtot)
      @param sstate (input) the state vector containing all values,
      scattered to processor 0. 
      @param time_data (input) a pointer to a struct that tipically is a time
      @return the double value corresonding to that dof */ 
  double get_dofval(const int & jeq,const double *sstate,
		     const TimeData *time_data) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a node/field  value from state vector.
      @author M. Storti
      @param node (input) the node number to get the row
      @param kdof (input) the field number to get the row
      @param sstate (input) the state vector
      @param ghost_vals (input) double array containing scattered ghost values
      @param time_data (input) a pointer to a struct that tipically is a time
      @param val (output) the double value corresonding to this
      node/field pair. */ 
  int get_nodal_value(const int & node,const int & kdof,double const *sstate,
		      double const *ghost_vals,const TimeData *time_data,double & val
		      ) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets a node/field  value from state vector.
      @author M. Storti
      @param node (input) the node number to get the row
      @param kdof (input) the field number to get the row
      @param sstate (input) the state vector (scattered to processor
      0. 
      @param time_data (input) a pointer to a struct that tipically is a time
      @param val (output) the double value corresonding to this
      node/field pair. */ 
  int get_nodal_value(const int & node,const int & kdof,
		      const double * sstate,const TimeData *time_data,
		      double & value) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets a fixation. 
      @author M. Storti
      @param node (input) the node number to set
      @param kdof (input) the field number to set
      @param val (input) the value to be set. */ 
  void set_fixation(int node,int kdof,double val);

  double value(const fixation_entry &fe, const TimeData *time_data) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns the unique number corresponding to the node/field
      pair. This is the inverse of nodf().
      @author M. Storti
      @param node (input) the node 
      @param kdof (input) the field 
      @return the unique number */ 
  int edof(const int node,const int field) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /**  Transforms nodf to (node,field). This is the inverse of edof().
      @author M. Storti
      @param edof (input) the nodf unique value
      @param node (output) the node 
      @param kdof (output) the field */ 
  void nodf(int edof,int &node,int &field) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns true if col j is null. 
      @author M. Storti
      @param j (input) the column index. */ 
  int col_is_null(int j);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets a constraint from the list of node/field/coefficients. 
      @author M. Storti
      @param constraint (input) list of node/field/coefficients. 
      @return 0/1 if the linear constraint has been detected to be
      linearly dependent with the preexisting constraints. */
  int set_constraint(const Constraint &constraint);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns the range of dofs that belong to a given processor. 
      The dofs in this processor (base 0) are dof1 <= dof <= dof2. 
      @author M. Storti
      @param myrank (input) the identifier of this processor. 
      @param dof1 (input) the start of the dof range. 
      @param dof2 (input) the end of the dof range. */ 
  void dof_range(int myrank,int &dof1,int &dof2);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Creates an MPI vector.
      @author M. Storti
      @param v (output) the vector to be created. */ 
  int create_MPI_vector(Vec &v);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Creates an MPI vector with the ghost values. 
      @author M. Storti
      @param v (output) the vector to be created. */ 
  int create_MPI_ghost_vector(Vec &v);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  Dofmap();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  ~Dofmap();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  int processor(int j) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  void freeze();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solves an equation of the form Q x = y for x. 
      Assumes that rank(Q) = n = dim(x). 
      @author M. Storti 
  */
  void solve(double *y,double *x);

  /** Makes #y=y+alpha*Q*x#. 
      @param y (input/output) the vector to add #Q*x# (size #nnod*ndof#). 
      @param x (input) the vector to multiply by #Q# (size #neq#).
      @param alpha (input) the scalar to scale the term #Q*x#. */ 
  void qxpy(double *x,double *y,double alpha);

  /** Makes #y=y+alpha*Q*x#. 
      @param y (input/output) the vector to add #Q'*x# (size #neq#). 
      @param x (input) the vector to multiply by #Q'# (size #nnod*ndof#).
      @param alpha (input) the scalar to scale the term #Q*x#. */ 
  void qtxpy(double *x,double *y,double alpha);

  int mult(Mat A,Vec x,Vec y);
};
#endif

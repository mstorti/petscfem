// -*- mode: c++ -*-
/*__INSERT_LICENSE__*/
//$Id: fastmat2.h,v 1.23 2002/12/04 03:14:12 mstorti Exp $

#ifndef FASTMAT2_H
#define FASTMAT2_H

#include <stdarg.h>

#include <utility>
#include <algorithm>
#include <vector>
#include <deque>
#include <cassert>

#include <newmatio.h>

#include "fastlib.h"
#include "readlist.h"

/// To be used in variable argument functions. Ex: fun(double,INT\_ARG\_LIST)
#define INT_ARG_LIST ARG_LIST(int,arg,0)
/// To access the elements in the INT\_ARG\_LIST.
#define READ_INT_ARG_LIST(indx) READ_ARG_LIST(arg,indx,0,EXIT)

/** Defines whether uses standard variadic argument lists or that
    defined rhtough readlist.h
*/
//#define USE_VAR_ARGS

#ifdef USE_VAR_ARGS
#define INT_VAR_ARGS ...
#else
#define INT_VAR_ARGS INT_ARG_LIST
#endif

#define FASTMAT_BASE_1

#ifdef FASTMAT_BASE_1
#define FASTMAT_BASE 1
#define LOCATION location1
#else
#define FASTMAT_BASE 0
#define LOCATION location0
#endif
 
//#define FASTMAT_DEBUG
//#define FM_DBG_COUNT

#ifdef FM_DBG_COUNT
#define FM_COUNT(count)   printf("FM count: %d\n",count)
#else
#define FM_COUNT(count)
#endif

#if 0
#ifdef FASTMAT_DEBUG
#define FMCHK if (ierr>0) printf("FastMat error at line %d, file: %s\n",__LINE__,__FILE__) 
#else
#define FMCHK 
#endif
#endif
  
/// Prints a FastMat2 matrix.
#define FMSHV(a) (a).print(#a)

//#define FM2_CACHE_DBG

//enum IndxOPT { LIST };

// The following is obsolete. Now it is defined in a `fastlib' vector
// as a FastVector template class.
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Fast version
#if 0
#define INDX_CHUNK_SIZE 10
class Indx {
public:
  inline Indx(void) {size_=0; storage=INDX_CHUNK_SIZE; 
  store=rigid_store; flexible_store=NULL;}; 
  inline Indx(const int m,const int n);
  inline ~Indx();
  void print(const char *s=NULL) const;
  int & operator[] (const int j) {return store[j];};
  const int & operator[] (const int j) const {return store[j];};
  inline int operator== (const Indx & indx) const;
  inline Indx & operator= (const Indx & indx);
  int size(void) const {return size_;};
  int push_back(const int j) {resize(size_+1); store[size_++]=j;};
  void reset() {size_=0;};
private:
  inline void resize(const int n);
  int rigid_store[INDX_CHUNK_SIZE];
  int *flexible_store;
  int *store;
  int storage;
  int size_;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Indx & Indx::operator= (const Indx & indx) {
  store = rigid_store;
  if (flexible_store) {
    delete[] flexible_store;
  }
  resize(indx.size_);
  size_ = indx.size_;
  for (int j=0; j<size_; j++) {
    store[j] = indx.store[j];
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline void Indx::resize(const int n) {
  if (n <= storage) {
    return;
  } else {
    storage += INDX_CHUNK_SIZE;
    if (storage<n) storage = n;
    int *new_flexible_store = new int[storage];
    for (int j=0; j<size_; j++) 
      new_flexible_store[j] = store[j];
    if (flexible_store) {
      delete[] flexible_store;
    }
    store = flexible_store = new_flexible_store;
  }
}
  
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline Indx::~Indx() {
  if(flexible_store) {
    delete[] flexible_store;
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline int Indx::operator== (const Indx & indx) const {
  if (size_!=indx.size_) return 0;
  for (int j=0; j<size_; j++) {
    if (store[j] != indx.store[j]) return 0;
  }
  return 1;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline Indx::Indx(const int m,const int n) {
  flexible_store=NULL;
  store=rigid_store;
  storage=INDX_CHUNK_SIZE;
  size_=0;
  resize(m);
  size_=m;
  for (int j=0; j<m; j++) {
    store[j]=n;
  }
}
#endif 

/// The type of indices used internally in FastMat2
typedef FastVector<int> Indx;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** An object of class perm stores the permutation of indices inside
    FastMat2 
*/
class Perm : public vector<int> {
public:
  Perm(void) {}; 
  Perm(const int m);
  void print(const char *s=NULL) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class IndexFilter : public Indx {
public:
  IndexFilter(const int dim_=0) : dim(dim_){};
  // adds a range to the filter
  void push(const int s,const int f=0,const int st=1);
  void print() const;
  int dim;
  // returns the absolute index corresponding to a given mapped index
  inline int abs_indx(const int j) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
typedef vector<IndexFilter> DimArray;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class LineCache {
public:
  vector<double *> linea,lineb;
  double *target,**starta,**startb;
  int linear,inca,incb;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
// Forward declarations
class FastMatCacheList;
typedef pair<FastMatCacheList *,int> FastMatCachePosition;

class FastMatSubCache {
public:
  virtual ~FastMatSubCache()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class FastMatCache {
public:
  FastMatCache();
  ~FastMatCache();
  // For one_to_one or in_place operations
  vector<double *> to_elems,from_elems;
  int nelems;
  double **pto,**pfrom,*to,*from;
  // For product applications (line,line)->double or line->double
  vector<LineCache> prod_cache;
  LineCache *line_cache_start;
  int nlines,line_size;
  vector<FastMatCacheList *> branch;
  Matrix *A,*B;
  FastMatSubCache *sc;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// This stores the counters for the various kind of operations
struct OperationCount {
  OperationCount() {get=0;put=0;mult=0;sum=0;div=0;abs=0;fun=0;};
  int get,put,mult,sum,div,abs,fun;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
//typedef vector<FastMatCache *> FastMatCacheList;
class FastMatCacheList : public vector<FastMatCache *> {
public:
  FastMatCacheList() {list_size=0;};
  // Operation count
  //  OperationCount op_count;
  //  void print_count_statistics();
  //  double operation_count(void);
  int list_size;
};

typedef double scalar_fun_with_args_t(double,void*);
typedef double scalar_fun_t(double);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Fast matrices. Supposed to be faster than Newmat
class FastMat2 {
public:
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Default constructor
      @author M. Storti
   */ 
  FastMat2(void);
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Constructor from indices (reading with varargs).
      @author M. Storti
      @param m (input) the number of indices. 
      @param INT\_VAR\_ARGS (input) the list of dimensions (m values)
  */ 
  FastMat2(const int m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Constructor from indices
      @author M. Storti
      @param dims\_ (input) the vector of dimensions.
  */ 
  FastMat2(const Indx & dims_);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Destructor
      @author M. Storti
  */ 
  ~FastMat2();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns the jd-th dimension
      @author M. Storti
      @param jd (input) returns the jd-th dimension
      @return the corresponding dimension
  */ 
  int dim(const int jd) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns the number of indices
      @author M. Storti
      @return the number of indices 
  */ 
  int n() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Flags if matrix was already defined
      @author M. Storti
      @param void (input) Used internally, flags if the matrix was
      already defined or not. 
      @return the correponding logical value
  */ 
  int is_defined(void) const {return defined;};

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Prints the matrix.
      @author M. Storti
      @param s (input) An optional string to print.
      @return a reference to the matrix.
   */ 
  void print(const char *s=NULL) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Prints dimensions of the matrix.
      @author M. Storti
      @param s (input) An optional string to print.
      @return a reference to the matrix.
  */ 
  void printd(char *s=NULL) ;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Add a range to the specified index filter. 
      @doc Adds the range (for j=start; j\<finish; j+=step) of indices
      to the current filter list. If start is undefined resets the
      filter for this index.
      @author M. Storti
      @param index (input) The index to which a range in the filter is
      added. 
      @param start (input) the start of the range
      @param finish (input) the end of the range. Def: start
      @param step (input) the step to be used. Def: 1
      @return a reference to the matrix.
  */ 
  FastMat2 & is(const int index,const int start=0,const int finish=0,
	    const int step=1);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Get filtered dims.
      @author M. Storti
      @param indx (output) The dimensions of the filtered matrix.
  */ 
  void get_dims(Indx & indx) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Adds a range to the row filter. 
      @doc Same as is() but acts on the first index. 
      @author M. Storti
      @param start (input) the start of the range
      @param finish (input) the end of the range
      @param step (input) the step to be used.
      @return a reference to the matrix.
  */ 
  FastMat2 & r(const int i,const int j=0,
	    const int step=1) 
    {return is(1,i,j,step);};


  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Adds a range to the column filter. Same as is() but acts on the
      second index.
      @author M. Storti
      @param start (input) the start of the range
      @param finish (input) the end of the range
      @param step (input) the step to be used.  
      @return a reference to the matrix.
  */
  FastMat2 & c(const int i,const int j=0,
	    const int step=1) 
    {return is(2,i,j,step);};

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets an index to a fixed value and reducing one dimension. 
      @doc The resulting filtered matrix is one dimension lower, and
      the corresponding index `indx' is set to `j'. 
      @author M. Storti
      @param indx (input) the index to be restricted.
      @param j (input) the value to which it is restricted. 
      @return a reference to the matrix.
   */ 
  FastMat2 & ir(const int indx,const int j=0);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets an index to mirror another, so that it creates a 
      mask that lets see the diagonal part of the matrix. 
      @author M. Storti
      @param j1 (input) the first index
      @param j2 (input) the second index
      @return a reference to the matrix.
  */ 
  FastMat2 & d(const int j1,const int j2);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Reset index filters for all dimensions. 
      @author M. Storti
      @return a reference to the matrix.
  */ 
  FastMat2 & rs();

  /** @name Operations on individual elements 
      @doc These operations act on individual elements of the matrix.
   */

  //@{
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets value at filtered position indx. 
      @author M. Storti
      @param indx (input) the vector of indices defining the position
      to be set. 
      @param val (input) the value to be set.
      @return a reference to the matrix.
  */ 
  FastMat2 & setel(const Indx & indx,const double val) {
    *location(indx) = val; };

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Adds value at filtered position indx. 
      @author M. Storti
      @param indx (input) the vector of indices defining the position
      to be set. 
      @param val (input) the value to be added.
      @return a reference to the matrix.
  */ 
  FastMat2 & addel(const Indx & indx,const double val) {
    *location(indx) += val; };

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets value at filtered position i,j,k... (INT\_VAR\_ARGS)
      @author M. Storti
      @param val (input) the value to be set.
      @param INT\_VAR\_ARGS (input) the indices of the position to be set.
      @return a reference to the matrix.
  */ 
  FastMat2 & setel(const double val,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Adds value at filtered position i,j,k... (INT\_VAR\_ARGS)
      @author M. Storti
      @param val (input) the value to be added.
      @param INT\_VAR\_ARGS (input) the indices of the position to be set.
      @return a reference to the matrix.
  */ 
  FastMat2 & addel(const double val,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Multiplies element at filtered position i,j,k... (INT\_VAR\_ARGS) by val. 
      @author M. Storti
      @param val (input) the coefficient to multiply.
      @param INT\_VAR\_ARGS (input) the indices of the position to be multiplied.
      @return a reference to the matrix.
  */ 
  FastMat2 & multel(const double val,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns the value at a given position.
      @author M. Storti
      @param i,j,k,l (input) the indices of the position
      @return the value at that position.
   */ 
  double get(INT_VAR_ARGS) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets value at filtered position 
      @author M. Storti
      @param indx (input) the vector of indices to retrieve element
      @return the value at that position.
   */ 
  double get(const Indx & indx) const {
    return *location(indx); };
  //@}

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Exchange indices (Obsolete)
      @author M. Storti
      @param i1 (input) One of the indices  to be exchanged 
      @param i2 (input) The other index to be exchanged 
      @return a reference to the matrix.
   */ 
  FastMat2 & exc(const int i1,const int i2);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Transpose (for matrix with two indices)
      @author M. Storti
      @return a reference to the matrix.
  */ 
  FastMat2 & t(void) {
    assert(dims.size()==2); return exc(1,2); };
 
  /** @name One to one operations.  
      These operations act on pair of
      matrices transforming from one to the other element by element.
  */
  //@{

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Copy matrix.
      @author M. Storti
      @param A (input) matrix to copy from
      @return a reference to the matrix.
  */ 
  FastMat2 & set(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Add matrix.
      @author M. Storti
      @param A (input) matrix to add
      @return a reference to the matrix.
  */ 
  FastMat2 & add(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Substract a matrix.
      @author M. Storti
      @param A (input) matrix to substract
      @return a reference to the matrix.
  */ 
  FastMat2 & rest(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** @memo Multiply (element by element) (like Matlab {\tt .*}). 
      @author M. Storti
      @param A (input) matrix to multiply
      @return a reference to the matrix.
  */ 
  FastMat2 & mult(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** @memo Divide matrix (element by element, like Matlab {\tt ./}). 
      @author M. Storti
      @param A (input) matrix to divide
      @return a reference to the matrix.
  */ 
  FastMat2 & div(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** @memo Reciprocal, #B.rcp(A,c)# is equivalent to Matlab (#B = c ./ A#). 
      @author M. Storti
      @param A (input) matrix to take the reciprocal.
      @return a reference to the matrix.
  */ 
  FastMat2 & rcp(const FastMat2 & A,double c=1.);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Axpy operation (element by element): {\tt (*this) += alpha * A}
      @author M. Storti
      @param A (input) matrix to add. 
      @param alpha (input) 
      @return a reference to the matrix.
  */ 
  FastMat2 & axpy(const FastMat2 & A,double alpha);
  //@}

  /** @name In-place operations. 
      @doc These operations perform an action on
      all the elements of a matrix.
  */
  //@{

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets all the element of a matrix to a constant value. 
      @author M. Storti
      @param val (input) the value to be set
      @return a reference to the matrix.
   */ 
  FastMat2 & set(const double val=0.);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Scale by a constant value.
      @author M. Storti
      @param val (input) the scale factor
      @return a reference to the matrix.
   */ 
  FastMat2 & scale(const double val);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Adds constant val
      @author M. Storti
      @param val (input) the value to be added
      @return a reference to the matrix.
  */ 
  FastMat2 & add(const double val);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Computes reciprocal of elements. In Matlab notation #A.rcp(c)# is
      equivalent to #A = c./A#.
      @author M. Storti
      @param c (input) scales the reciprocal
      @return a reference to the matrix.
  */ 
  FastMat2 & rcp(const double c=1.);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Apply a function to all elements
      @author M. Storti
      @param function (input) the function to be applied
      @return a reference to the matrix.
  */ 
  FastMat2 & fun(scalar_fun_t *function);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Apply a function with optional arguments to all elements
      @author M. Storti
      @param function (input) the function to be applied
      @param user_args (input) additional arguments to the scalar function
      @return a reference to the matrix.
  */ 
  FastMat2 & fun(scalar_fun_with_args_t *function,void *user_args);
  //@}

  /** @name Miscellaneous utilities. */
  //@{
  /** Cross product of `a' and `b'. (all must be vectors of length 3).
      @author M. Storti
      @param a (input) first vector
      @param b (input) second vector
      @return a reference to the matrix.
  */ 
  FastMat2 & cross(const FastMat2 & a,const FastMat2 & b);
  /** Sets to the Levi-Civita density tensor. This is a third 
      order tensor $\epsilon_{ijk}$ ($3\times3\times3$) 
      such that $\epsilon_{ijk}$ is 1 if $ijk$ is 
      an even permutation of 123, -1 if it is an odd permutation 
      and 0 otherwise. The vector cross-product and determinant p
      can be computed in terms of this tensor. 
      @author M. Storti
      @param a (input) the value to be set in the diagonal
      @return a reference to the matrix.
  */ 
  FastMat2 & eps_LC();

  /** Sets to a multiple of the identity matrix.
      @author M. Storti
      @param a (input) the value to be set in the diagonal
      @return a reference to the matrix.
  */ 
  FastMat2 & eye(const double a=1.);
  //@}

  /** @name Product and contraction operation. */
  //@{
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Contraction operation. Contracts pairs of indices. (generalized trace). 
      @author M. Storti
      @param A (input) The matrix to take the trace
      @param i,j,k,l.... (input) The indices to contract and remap indices.
      @return a reference to the matrix.
  */ 
  FastMat2 & ctr(const FastMat2 & A,const int m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Product operation. {\tt (*this) = A * B} (generalized by index contraction)
      @author M. Storti
      @param A (input) first natrix
      @param B (input) second natrix
      @param i,j,k,l... (input) indices to contract and remap.
      @return a reference to the matrix.
  */ 
  FastMat2 & prod(const FastMat2 & A,const FastMat2 & B,const int
		  m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Kronecker product (also called Schur product)
      @doc If A is n x m and B is p x q, then kron returns a matrix
      which is np x mq, and where each p x q block is proportional to
      B. Only implemented for two-dimensional matrices so far.
      @author M. Storti
      @param A (input) first matrix 
      @param B (input) second matrix 
      @return a reference to the matrix.
   */ 
  FastMat2 & kron(const FastMat2 & A,const FastMat2 & B);
  //@}

  /** @name Sum operations (sum over indices), #max#, #sum_#, ...
      These operations give a reduced matrix by applying some reduce
      operation (sum, max, min etc...) over all the contracted
      indices. For instance if A has 4 indices then {\tt B.sum(A,2,-1,-1,1)}
      is equivalent to $B_{ij} = \sum_{kl} A_{jkli}$
   */
  //@{

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & sum(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum of squares over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & sum_square(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum of absolute values all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & sum_abs(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Norm p of matrix (per column). 
      $(\sum_j |a_j|^p)^{\frac 1p}$
      @author M. Storti
      @param A (input) matrix to contract
      @param p (input) exponent of norm
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & norm_p(const FastMat2 & A,double p,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Norm #p# of matrix (per column) for #p# integer (more efficient
      than the version for #p# real)
      $(\sum_j |a_j|^p)^{\frac 1p}$
      @author M. Storti
      @param A (input) matrix to contract
      @param p (input) exponent of norm
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & norm_p(const FastMat2 & A,int p,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Minimum over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & min(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Maximum over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & max(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Min of absolute value  over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & min_abs(const FastMat2 & A,const int m=0,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Min of absolute value over all selected indices.
      @author M. Storti
      @param A (input) matrix to contract
      @param i,j,k,l... (input) indices that define indices to be contracted
      @return a reference to the matrix.
   */ 
  FastMat2 & max_abs(const FastMat2 & A,const int m=0,INT_VAR_ARGS);
  //@}


  /** @name Sum operations over all indices, max, sum\_, ...
      These operations give a scalar by applying some reduced
      operation (sum, max, min etc...) over all 
      indices. For instance if A has 4 indices then {\tt b = sum\_all(A)} returns
      a double $b = \sum_{kl} A_{jkli}$
   */
  //@{

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double sum_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum of squares over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double sum_square_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sum of absolute values over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double sum_abs_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Norm p over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double norm_p_all(const double p) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Norm p over all indices for integer p.
      @author M. Storti
      @return the result of the operation
  */ 
  double norm_p_all(const int p) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Minimum  over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double min_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Maximum over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double max_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Minimum absolute value over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double min_abs_all() const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Maximum absolute value over all indices.
      @author M. Storti
      @return the result of the operation
  */ 
  double max_abs_all() const;
  //@}

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Reshapes a matrix.
      New and old dimensions must be compatible. 
      @author M. Storti
      @param i,j,k,l (input) the new dimensions
      @return a reference to the matrix.
  */ 
  FastMat2 & reshape(const int m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Resizes the matrix (destroying the information).
      New and old dimensions might not be the same. 
      @author M. Storti
      @param i,j,k,l (input) the new dimensions
      @return a reference to the matrix.
   */ 
  FastMat2 & resize(const int m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Resizes to one dimension ans zero elements.
   */ 
  FastMat2 & clear();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Sets vector to diagonal part 
      @author M. Storti
      @param A (input) the matrix to take the diagonal
      @param i,j,k,l... (input) the indices defining the operation
      @return a reference to the matrix.
  */ 
  FastMat2 & diag(FastMat2 & A,const int m,INT_VAR_ARGS);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** @name Export/Import operations
   */ 
  //@{
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Copies to argument from Newmat matrix.
      @author M. Storti
      @param A (input) the matrix to be copied
      @return a reference to the matrix.
  */ 
  FastMat2 & set(const Matrix & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Copy from array of doubles.
      @author M. Storti
      @param a (input) the array of indices from where to copy
      @return a reference to the matrix.
  */ 
  FastMat2 & set(const double *a);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** exports to a double vector
      @author M. Storti
      @param a (output) array doubles to where export
      @return a reference to the matrix.
  */ 
  const FastMat2 & export_vals(double *a) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** exports to a double vector
      @author M. Storti
      @param a (output) array doubles to where export
      @return a reference to the matrix.
  */ 
  FastMat2 & export_vals(double *a) ;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Exports to a Newmat matrix.
      @author M. Storti
      @param A (output) the Newmat matrix to where export to
      @return a reference to the matrix.
  */ 
  const FastMat2 & export_vals(Matrix & A) const;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Exports to a Newmat matrix.
      @author M. Storti
      @param A (output) the Newmat matrix to where export to
      @return a reference to the matrix.
  */ 
  FastMat2 & export_vals(Matrix & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Returns a pointer to the start of the storage matrix. 
      @author M. Storti
      @returns as mentioned
  */ 
  double *storage_begin();
  //@}
  
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Determinant
      @author M. Storti
      @return the determinant of the matrix.
   */ 
  double det(void) const;

  /** For a matrix of size #A# of #(n-1)* n# computes the determinant
      #sqrt(det(A*A'))#. This is useful when integrating on surfaces
      in 3D and lines in 2D. If #nor# is not null then computes #nor#
      as the normal to the surface. 
      @param nor (ouput) the normal to the surface
      @return determinant of the surface Jacobian
  */ 
  double detsur(FastMat2 *nor=NULL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** The inverse of a matrix. {\tt (*this) = inverse of A}
      @author M. Storti
      @param A (input) the matrix to take the inverse 
      @return a reference to the matrix.
  */ 
  FastMat2 & inv(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solve the eigenvalue problem for non-symmetric matrix #A# and 
      compute right eigenvectors. (see 
      #eig(const FastMat2 & A, FastMat2 *VR=NULL,...)# below).
      @author M. Storti
      @param A (input) the matrix to take the eigenvalues
      @param VR (output) the matrix with the right eigenvectors
      @returns a reference to a vector containing the eigenvalues. 
  */ 
  FastMat2 & eig(const FastMat2 & A, FastMat2 &VR);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solve the eigenvalue problem for non-symmetric matrix #A# and 
      compute right and left eigenvectors. (see 
      #eig(const FastMat2 & A, FastMat2 *VR=NULL,...)# below).
      @author M. Storti
      @param A (input) the matrix to take the eigenvalues
      @param VR (output) the matrix with the right eigenvectors
      @param VL (output) the matrix with the left eigenvectors
      @returns a reference to a vector containing the eigenvalues. 
  */ 
  FastMat2 & eig(const FastMat2 & A, FastMat2 &VR, FastMat2 &VL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solve the eigenvalue problem for non-symmetric matrix #A#.
      #A# has to be square. If #A# is of size #m x m# then it returns
      the eigenvalues in a matrix of size #2 x m# containing the real and
      imaginary part of the eigenvalues. Both right and left eigenvectors
      can be computed independently. Based on the #DGEEV# #LAPACK# routine. 
      @author M. Storti
      @param A (input) the matrix to take the eigenvalues
      @param VR (output) the matrix with the right eigenvectors
      @param VL (output) the matrix with the left eigenvectors
      @param crev (input) compute the right eigenvectors flag
      @param clev (input) compute the left eigenvectors flag
      @returns a reference to a vector containing the eigenvalues. 
  */ 
  FastMat2 & eig(const FastMat2 & A, FastMat2 *VR=NULL, FastMat2 *VL=NULL, 
		 int crev=0, int clev=0);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solve the eigenvalue problem for symmetric matrix #A#.
      Note: #A# is NOT checked for symmetry. #A# has to be square. 
      This routine may be more faster and accurate tha the 
      non-symmetric version. Based on the #DSYEV# #LAPACK# routine. 
      @author M. Storti
      @param A (input) the matrix to take the eigenvalues
      @returns a reference to a vector containing the eigenvalues. 
  */ 
  FastMat2 & seig(const FastMat2 & A);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Solve the eigenvalue problem for symmetric matrix #A#.
      Note: #A# is NOT checked for symmetry. #A# has to be square. 
      Based on the #DSYEV# #LAPACK# routine. 
      @author M. Storti
      @param A (input) the matrix to take the eigenvalues
      @param V (output) the matrix with the eigenvectors
      @param compute_eigen_vectors (input) flags whether
      to compute the eigenvectors or not
      @returns a reference to a vector containing the eigenvalues. 
  */ 
  FastMat2 & seig(const FastMat2 & A, FastMat2 &V,int compute_eigen_vectors=1);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Converts to double, for zero dimension matrices.
      @author M. Storti
      @return a double that is the only element of the matrices.
   */ 
  operator double() const;

  // Converts to newmat
  // friend void FM2_2_NM(Matrix & M,const FastMat2 & M_FM2);

  // Sets current cache list
  //  static void set_cache_list(FastMatCacheList & cache_list_);

  /** @name Static cache operations
   */
  //@{

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Activates use of the cache
      @author M. Storti
      @param cache_list_ (input) the cache list root to activate.
  */ 
  static void activate_cache(FastMatCacheList *cache_list_=NULL);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Deactivates use of the cache
      @author M. Storti
  */ 
  static inline void deactivate_cache(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Resets the cache.
      To be used after each iteration loop.
      @author M. Storti
  */ 
  static void reset_cache(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Voids the cache.
      Frees the memory used by the cahe list after processing.
      @author M. Storti
  */ 
  static void void_cache(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Creates a branch point.
      @author M. Storti
      @param (input)
  */ 
  static void branch(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Follows a branch. 
      @author M. Storti
      @param j (input) the number of branch to follow.
  */ 
  static void choose(const int j);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Leaves the current branch.
      @author M. Storti
  */ 
  static void leave(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Computes the total number of operations in the cache list.
      Currently all operations counts as one. In the future we will
      give weights to each type of operation. 
      @author M. Storti
      @return the total number of operations.
  */ 
  static double operation_count(void);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Print statistics about the number of operations of each type in
      the current cache list.
      @author M. Storti
  */ 
  static void print_count_statistics();

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Gets actual position in cache. 
      @author M. Storti
  */ 
  static void get_cache_position(FastMatCachePosition & pos);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Jumps to a given position in the cache-list.
      @author M. Storti
  */ 
  static void jump_to(FastMatCachePosition &pos);

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /** Recomputes the was\_cached variable
      @author M. Storti
  */ 
  static void resync_was_cached(void);
  //@}

private:

  /// Total storage. Should be the product of `dims'. 
  int storage;
  /// The array of values. Should be of size `storage'. 
  double *store;
  /// dimensions
  DimArray dims;
  /// Pointer to dims, as an array
  IndexFilter *dims_p;
  /// Size of dims array
  int n_dims;
  /// Permutation of indices(transposing)
  Perm perm;
  /// Fixed indices 
  Indx set_indx;
  /// auxiliary vector to put absolute indices
  Indx absindx;
  /// was the matriz defined? 
  int defined;
  /// Root of cache lists
  static FastMatCacheList *cache_list_root;
  /// Current list of caches
  static FastMatCacheList *cache_list;
  /// Cache-list stack
  static vector<FastMatCachePosition> cache_list_stack;
  /// Position in cache\_list
  static int position_in_cache;
  /// begin of cache\_list
  static FastMatCache **cache_list_begin;
  /// size of cache list
  static int cache_list_size;
  /// Use cache?
  static int use_cache;
  /// Was computed this cache list
  static int was_cached;
  /// save was\_cached if use deactivate\_cache()
  static int was_cached_save;
  /// Operation count
  static OperationCount op_count;
  int comp_storage_size(const Indx & indx) const;
  /// creates storage and freezes dimensions
  void define_matrix(void);
  /// returns value at position i,j
  double val(const int i,const int j) const;
  /// returns address of absolute  position indx[0],indx[1],...
  double *location_abs(const Indx & indx) const;
  /// returns address of filtered position indx[0],indx[1],...
  double *location(const Indx & indx) const;
  // gets value at absolute position indx 
  void set_abs(const int i, const int j,const double val);
  /// used in constructors
  void create_from_indx(const Indx & dims_);
  /// auxiliary.  prints matrices with 2 indices.
  void print2(const Indx & indxp,const Indx & fdims) const;
  /// auxiliary.  prints matrices with 1 indices.
  void print1(const Indx & indxp,const Indx & fdims) const;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline void FastMat2::deactivate_cache(void) {
  was_cached_save=was_cached;
  use_cache=0; 
  was_cached=0; 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// For 2 indices matrices
class FMatrix : public FastMat2 {
public:
  FMatrix() : FastMat2() {};
  FMatrix(const int m, const int n) : FastMat2(2,m,n) {};
  FMatrix(const int m) : FastMat2(1,m) {};
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// increment an index
int inc(Indx & indx,const Indx &dims);

/// Read a list of integers from variable arguments in a vector 
void read_int_list(const int m,va_list v,Indx *indx);

/// prints a FastMat2 matrix for debugging.
#define FSHV(name) (name).print(#name ": ")

/// Obtains the amount of memory used by a matrix from its vector of
/// indices
int mem_size(const Indx & indx);

double int_pow(double base,int exp);
#endif

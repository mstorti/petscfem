// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse.h,v 1.26 2001/11/09 03:05:42 mstorti Exp $
#ifndef SPARSE_H
#define SPARSE_H

#include <cmath>
#include <cstdio>

#include <map>
#include <vector>
#include <algorithm>

#include <SRC/util.h>
#include <SRC/dsp_defs.h>

#include <sles.h>

#include "randomg.h"

using namespace Random;

namespace Sparse {

  const double MaxDouble=1.e300;
  const double MinDouble=-1.e300;

  class ScalarFunObj {
  public:
    virtual ~ScalarFunObj() =0;
    virtual double fun(double v) const=0;
  };

  class Scale : public ScalarFunObj {
  public:
    double c;
    ~Scale() {};
    double fun(double v) const {return c*v;};
  };

  extern Scale scale_fun_obj;

  typedef double ScalarFunD(double v,void * user_data = NULL);
  typedef double ScalarFun(double v);

  class ScalarFunWrapper : public ScalarFunObj {
  public:
    ~ScalarFunWrapper() {};
    void * user_data;
    ScalarFunD *sfd;
    ScalarFun *sf;
    virtual double fun(double v) const;
    ScalarFunWrapper() : sfd(NULL), sf(NULL) {};
  };

  extern ScalarFunWrapper scalar_fun_wrapper;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class BinAssoc {
  public:
    double null;
    virtual ~BinAssoc() =0;
    // `v' is the "cumulated" value!! This may help in improve
    // efficiency. 
    virtual double op(double v,double w) const = 0;
  };

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  class Sum : public BinAssoc {
  public:
    ~Sum() {}
    double op(double v,double w) const {return v+w;};
    Sum() {null=0;};
  };

  extern Sum sum_bin_assoc;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class SumAbs : public BinAssoc {
  public:
    ~SumAbs() {}
    double op(double v,double w) const {return v + fabs(w);};
    SumAbs() {null=0;};
  };

  extern SumAbs sum_abs_bin_assoc;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class Max : public BinAssoc {
  public:
    ~Max() {};
    double op(double v,double w) const {return (v>w? v : w);};
    Max() {null=MinDouble;};
  };

  extern Max max_bin_assoc;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class MaxAbs : public BinAssoc {
  public:
    ~MaxAbs() {};
    double op(double v,double w) const {
      double a = fabs(w); return (v>a? v : a);};
    MaxAbs() {null=0.;};
  };

  extern MaxAbs max_abs_bin_assoc;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class Min : public BinAssoc {
  public:
    ~Min() {};
    double op(double v,double w) const {return (v<w? v : w);};
    Min() {null=MaxDouble;};
  };

  extern Min min_bin_assoc;

  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 

  class Accumulator {
  public:
    virtual ~Accumulator() =0;
    virtual void accum(double &v,double w)=0;
  };

  class SumSq : public Accumulator {
  public:
    ~SumSq() {};
    void accum(double &v,double w) {v += w*w;};
  };

  extern SumSq sum_sq_accum;

  class SumPow : public Accumulator {
  public:
    double n;
    ~SumPow() {};
    void accum(double &v,double w) {v += pow(fabs(w),n);};
  };

  extern SumPow sum_pow_accum;

  class Mat;
      
  class GenVec {
  protected:
    /// Value of those elements that are not represented
    static double not_represented_val;
  public:
    /// Pure virtual class
    virtual ~GenVec()=0;
    /// Return length of the vector 
    virtual int length() const =0;
    /// Get element at specified position
    virtual double get(int j) const =0;
    /// Set element at position j
    virtual GenVec & set(int j,double v)=0;
    /// Matrix vector product
    virtual GenVec & prod(const Mat & a,const GenVec & v);
    /// print elements
    virtual void print(const char *s = NULL)=0;
    /// print elements
    virtual void print_f(const char *s = NULL);
    /// export to array
    virtual void get(double *s) const;
    /// Resize vector
    virtual GenVec & resize(int n)=0;
    /// export internal array
    virtual const double * get() const {assert(0);};
    /// generic copy vector
    virtual GenVec & set(const GenVec &v);
    
  };

  class FullVec : public vector<double>, public GenVec {
  public:
    /// Constructor
    FullVec() {};
    /// Constructor from GenVec
    FullVec(GenVec &);
    /// Copy from GenVec
    //    FullVec & set(GenVec & v);
    /// Return length of the vector 
    int length() const {return size();};
    /// Get element at specified position
    double get(int j) const {return (*this)[j];};
    /// Set element at position j
    FullVec & set(int j,double v) {(*this)[j] = v; return *this;};
    /// print elements
    void print(const char *s = NULL) {print_f(s);};
    /// export internal array
    const double * get() const {return begin();};
    /// Resize vector
    GenVec & resize(int n) {vector<double>::resize(n); return *this;};
  };

  class Indx : public vector<int> {
  public:
    Indx(int n,int m,int k=1);
  };
  typedef Indx::iterator IndxIt;
  typedef Indx::const_iterator IndxCIt;

  typedef map<int,double>::iterator VecIt;
  typedef map<int,double>::const_iterator VecCIt;
  typedef pair<int,double> VecP;

  /// A simple sparse vector class (0 indexed). 
  class Vec : public map<int,double>, public GenVec {

    /// Length of the sparse vector
    int len;
    /// Flag indicating where you can add values past the specified length or not 
    int grow_m; 
    /// Get element at specified position (do not check bounds)
    double get_nc(int j) const;
      
  public:

    friend class Mat;
    /// Constructor from the length
    Vec(int l=0) : grow_m(1) {len=l;};
    /// Return length of the vector 
    int length() const {return len;};
    /// Constructor from another vector
    Vec(const Vec &v) {*this = v;};

    /* Insert contents of vector v at position I.
	Length of vector is that of the maximum index in I. 
    */
    // Vec(const Indx &I,const Vec &v);
    /// Get element at specified position
    double get(int j) const;
    /// Get a subvector of elements at position I
    void get(const Indx &I,Vec &v) const; 

    /// Copy
    Vec & copy(const Vec &v) {*this = v; return *this;}; 

    /// Set element at position j
    Vec & set(int j,double v);
    /// Set elements at subvector at position I. w[I] = v
    Vec & set(const Indx & I,const Vec & v);
    /// Set vector to elements of v at position I. w = v[I]
    Vec & set(const Vec & v,const Indx & I);
    /// Set vector elements at I to elements f of v at position J. w[I] = v[J]
    Vec & set(const Indx & I,const Vec & v,const Indx & K);
    /// Set to row of a matrix
    Vec & set(const Mat & a,int j);
    /// Set vector to k-th column of matrix a. w = a(:,k)
    Vec & setc(const Mat &a,int k);

    /// Scale elements 
    Vec & scale(double c);

    /// Sets w += a * v
    Vec & axpy(double a,const Vec & v);
    /// Sets w[I] = c * w[I] + a * v
    Vec & axpy(double c,const Indx & I,double a,const Vec & v);

    /// Dot product
    double dot(const Vec & w) const;
    /// Dot product with generic vector
    double dot(const GenVec & w) const;

    /// print elements (generic version)
    void print_g(int l,const char * s,const char * psep, const char * isep,
		 const char * lsep) const;
    /// print elements
    void print(const char *s = NULL);
    /// Set mode if can grow automatically or not
    Vec & grow(int g) { grow_m=g; return *this;};
    /// Clears all elements
    Vec & clear() {map<int,double>::clear(); return *this;}
    /// Resize vectors, truncates elements if greater than this value
    Vec & resize(int n);
    /// Flags if the vector is empty or not
    int empty() const;
    /// Purge elements below a given tolerance value
    Vec & purge(double tol = 1e-10);
    /// Fill with random values
    Vec & random_fill(double fill=0.1,Generator & g = uniform);
    /// Apply a function object to all elements
    Vec & apply(const ScalarFunObj & fun);
    /// Apply a scalar function with args to all elements
    Vec & apply(ScalarFunD *fun,void * user_data = NULL);
    /// Apply a scalar function with no args to all elements
    Vec & apply(ScalarFun *fun);

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    /// perform an asscociative function on all non-null elements
    double assoc(BinAssoc & op) const;

    /// Sum of all elements
    double sum() const {return assoc(sum_bin_assoc);} ;
    /// Sum of absolute value of all elements
    double sum_abs() const {return assoc(sum_abs_bin_assoc);} ;
#if 0
    /// Sum of squares of all elements
    double sum_sq() const {return assoc(sum_sq_bin_assoc);} ;
    /// Sum of power of all absolute value of all elements
    double sum_pow(double n) const {sum_pow_bin_assoc.n = n;
    return assoc(sum_pow_bin_assoc);} ;
#endif
    /// Max of all elements
    double max() const {return assoc(max_bin_assoc);} ;
    /// Max of absolute value of all elements
    double max_abs() const {return assoc(max_abs_bin_assoc);} ;
    /// Min of all elements
    double min() const {return assoc(min_bin_assoc);} ;

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    /// perform an asscociative function on all non-null elements
    void accum(double &v,Accumulator & acc) const;

    /// Sum of squares of all elements
    double sum_sq() const {double v=0; accum(v,sum_sq_accum); return v;} ;
    /// Sum of power of all absolute value of all elements
    double sum_pow(double n) const {sum_pow_accum.n = n;
    double v; accum(v,sum_pow_accum); return v;};
  };
  
  typedef map<int,Vec>::iterator RowIt;
  typedef map<int,Vec>::const_iterator RowCIt;
  typedef pair<int,Vec> RowP;
  typedef pair<const int,Vec> RowCP;

  // The Sparse::Mat class is implemented as a finite state machine in
  // the sense proposed by Robert Martin. A base class, in this case
  // MatFSMContext implements all basic actions. SMC implements a
  // series of derived classes (for instance MatFSMfilledState), one
  // for each state of the machine, which contains the logic of the
  // machine. The user callable class, contains static members for
  // each state and a pointer to the static member that corresponds to
  // the actual state of the machine. The state of the machine
  // corresponds to what static member points this pointer. 

  // I successfully added some more complexity to this :-) . The user
  // callable class Mat contains an instance of the finite state
  // machine FSM class MatFSM.  Each method of Mat calls (or not) to a
  // FSM event. When this is traduced to an action to the base
  // MatFSMContext class this is executed on the Mat class via a
  // pointer `matrix_p' contained in the Context class. This
  // simplifies the fact in order that a event can trigger other
  // events. (In the proposal of R. Martin this is done otherwise). 

  // Actions appear in the definition of the FSM context class
  // header and definition. The definition is simply to call the
  // action on the `matrix_p' pointer .

#define FSM_ACTION_DECL(action) void action()

#if 1

#define FSM_ACTION_DEF(action) void MatFSMContext::action()	\
  {matrix_p->action();}

#else
  // For debugging
#define FSM_ACTION_DEF(action)			\
void MatFSMContext::action() {			\
    matrix_p->action();				\
    printf("IN ACTION: " #action		\
	   ", I'M %p\n",this);			\
}

#endif

#define FSM_ACTIONS				\
  FSM_OP(clear);				\
  FSM_OP(clean_mat);				\
  FSM_OP(clean_factor);				\
  FSM_OP(fact_and_solve);			\
  FSM_OP(solve_only)

  class MatFSMContext {
  public:
    Mat * matrix_p;
    MatFSMContext() {};
#undef FSM_OP
#define FSM_OP(action) FSM_ACTION_DECL(action)
    FSM_ACTIONS;
    void FSMError(const char *e,const char *s) { 
      printf("Not valid \"%s\" event in state \"%s\"",e,s);
    }
  };

#include "matFSM.h"

    // Simple sparse matrix class. 
  class Mat : public map< int, Vec >  {

    /// Dimensions
    int nrows,ncols;
    /// Flag indicating where you can add values past the specified dimensions 
    int grow_m; 
    /// Value of those elements that are not represented
    static double not_represented_val;

    void init_fsm(Mat *) {fsm.matrix_p = this;};
  public:

    static Mat *dispatch(char *opt);

    double *b;

    friend class MatFSM;
    MatFSM fsm;

    friend class GenVec;
    friend class Vec;

    /// Constructor from the length
    Mat(int m=0,int n=0) : grow_m(1), nrows(m), ncols(n)
			   { init_fsm(this);};

    /// Destructor (invokes clear() FSM)
    ~Mat() {clear();}

    /// Return row dimension
    int rows() const {return nrows;};
    /// Return column dimension
    int cols() const {return ncols;};
    /// Constructor from another vector
    Mat(const Mat &a) {*this = a; init_fsm(this);};

    /// Get element at specified position: v = w(j,k)
    double get(int j,int k) const;
    /// Get row: v = w(j,:)
    void getr(int j,Vec & v) const;
    /// Set element at position j: w(j,k) = v
    Mat & set(int j,int k,double v);

    /// Set row j: w(j,:) = v;
    Mat & setr(int j,Vec &v);
    /// set to J rows of a: w = a(J,:)
    Mat & setr(const Mat & a,Indx &J);

    /// Set col j: w(:,j) = v;
    Mat & setc(int j,const Vec &v);
    /// set to J cols of a: w = a(:,J)
    Mat & setc(const Mat & a,Indx &J);

    /// Set to multiple of Identity matrix 
    Mat & id(double a=1.);
    /// Diagonal from a vector
    Mat & diag(const Vec & v);
    /// Fill with random values
    Mat & random_fill(double fill=0.1,Generator & g=uniform);
    /// Product of sparse matrices
    Mat & prod(const Mat & a, const Mat & b,double c=1.);

    /// Apply a function object to all elements
    Mat & apply(const ScalarFunObj & fun);
    /// Apply a scalar function with args to all elements
    Mat & apply(ScalarFunD *fun,void * user_data = NULL);
    /// Apply a scalar function with no args to all elements
    Mat & apply(ScalarFun *fun);

    /// Sets w += a * v
    Mat & axpy(double c,const Mat & a);

    /// Scale elements 
    Mat & scale(double c) {scale_fun_obj.c = c; return apply(scale_fun_obj);};

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    /// perform an asscociative function on all non-null elements
    double assoc(BinAssoc & op) const;

    /// Sum of all elements
    double sum() const {return assoc(sum_bin_assoc);} ;
    /// Sum of absolute value of all elements
    double sum_abs() const {return assoc(sum_abs_bin_assoc);} ;
    /// Max of all elements
    double max() const {return assoc(max_bin_assoc);} ;
    /// Max of absolute value of all elements
    double max_abs() const {return assoc(max_abs_bin_assoc);} ;
    /// Min of all elements
    double min() const {return assoc(min_bin_assoc);} ;

    //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
    /// perform an asscociative function on all non-null elements
    void accum(double &v,Accumulator & acc) const;

    /// Sum of squares of all elements
    double sum_sq() const {double v=0.; accum(v,sum_sq_accum); return v;} ;
    /// Sum of power of all absolute value of all elements
    double sum_pow(double n) const {double v=0.; sum_pow_accum.n = n;
    accum(v,sum_pow_accum); return v;};

    /// print elements (sparse version)
    void print(const char *s = NULL) const;
    /// print elements (full version)
    void print_f(const char *s = NULL) const;

    /// Resize vectors, truncates elements if greater than this value
    Mat & resize(int m,int n);
    /// Clears all elements
    Mat & clear();
    /// Set mode if can grow automatically or not
    Mat & grow(int g) { grow_m=g; return *this;};
    /// Flags if the vector is empty or not
    int empty() const;
    /// Number of non null elements
    int size() const;
    /// Solve the linear system 
    void solve(FullVec &b);
    /// Solve the linear system 
    void solve(double *b);

    /// FSM actions
//  #undef FSM_OP
//  #define FSM_OP(action) FSM_ACTION_DECL(action)
//      FSM_ACTIONS;
    virtual void clean_factor()=0;
    void clean_mat();
    virtual void solve_only()=0;
    virtual void fact_and_solve()=0;

  };

}

#endif

// -*- mode: c++ -*-
//__INSERT_LICENSE__
//$Id: fastmat.h,v 1.2 2001/04/01 01:35:06 mstorti Exp $

#ifndef FASTMAT_H
#define FASTMAT_H

#include <cassert>
#include "vecmacros.h"

#define FASTMAT_BASE_1

#ifdef FASTMAT_BASE_1
#define FASTMAT_BASE 1
#define LOCATION location1
#else
#define FASTMAT_BASE 0
#define LOCATION location0
#endif
 
#define FASTMAT_DEBUG
//#define FM_DBG_COUNT

#ifdef FM_DBG_COUNT
#define FM_COUNT(count)   printf("FM count: %d\n",count)
#else
#define FM_COUNT(count)
#endif

#ifdef FASTMAT_DEBUG
#define FMCHK if (ierr>0) printf("FastMat error at line %d, file: %s\n",__LINE__,__FILE__) 
#else
#define FMCHK 
#endif
  
#define FMSHV(a) (a).print(#a)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/// Fast matrices. 
class FastMat {
public:
  /** Default constructor. Constructor from dimensions and from 
      optional vector
  */
  FastMat(int m=0,int n=0,double *s = NULL);
#if 0
  /// Constructor from a Newmat matrix
  FastMat(const Matrix & B);
#endif
  /// destructor
  ~FastMat();
  /// set to the identity matrix
  void eye(int m);
  /// sets correct size (only if not set previously)
  void set_size(int m,int n);
  /// sets element i,j to val
  void set(int i,int j,double val);
  /// sets submatrix  i1 to i2, j1 to j2 to val
  void set(int i1,int i2, int j1, int j2,double *val);
  /// sets submatrix  i1 to i2, j1 to j2 to matrix B
  void set(int i1,int i2, int j1, int j2,const FastMat & B);
  /// adds val to element i,j
  void add(const int i,const int j,const double val);
  /// adds submatrix in array val to submatrix i1 to i2, j1 to j1
  void add(int i1,int i2, int j1, int j2,double *val);
  /** adds submatrix from another matrix in array
      val to submatrix i1 to i2, j1 to j1
  */
  void add(int i1,int i2, int j1, int j2,const FastMat & B);
  /// copies from array. Only n first elements if n\\neq 0
  void set(double *val,int n=0);
  /// copies from other matrix. Only n first elements if n\\neq 0
  void set(const FastMat & A,int n=0);
  /// sets to all elements to a constant 
  void set(double val);
  /// returns element i,j (obsolete!!)
  void get(const int i,const int j,double *val) const;
  /// returns element i,j
  void get(const int i,const int j,double & val) const;
  /// returns element i,j
  double get(const int i,const int j) const {
    assert(defined && i-FASTMAT_BASE<m && j-FASTMAT_BASE<n );
    return *(LOCATION(i,j));
  };
  /// returns block (i1,i2,j1,j2) of matrix A a *this
  void get(const int i1,const int i2,const int j1,const int j2,const FastMat &A);
  /// prints the matrix
  void print(char *s);
  /// copies from matrix to external array. only first n elements if n.neq.0
  void copy(double *val,int n=0) const;
  /// copies a row to *this
  void row(const FastMat & B,int i);
  /// copies some rows to *this
  void rows(const FastMat & B,int i1,int i2);
  /// copies some columns to *this
  void column(const FastMat & B,int j);
  /// copies some columns to *this
  void columns(const FastMat & B,int j1,int j2);
  /// return pointer to location (i,j) (base 1)
  double *location1(const int i,const int j) const { return store+(i-1)*n+j-1;};
  /// return pointer to location (i,j) (base 0)
  double *location0(const int i,const int j) const { return store+i*n+j;};
  /// reshapes the matrix
  void reshape(int m, int n);
  /// converts to scalar
  int as_scalar(double & val) const {assert(m==1 || n==1); val=*store;}
  /// sum of all terms
  double sum() const;
  /// sum of squares
  int sum_square(double & val) const;
  /// transpose
  int transpose(const FastMat & A);
  /// scales the matrix
  int scale(const double c);
  /// Trace of product
  int trace_of_product(const FastMat & B,double & trace) const;
  /// Norm 1 
  double norm1() const;
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  //  TO BE SET AS PRIVATE
  //---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
  /// the storage
  double *store;
  /// matrix dimensions
  int m,n;
  /// number of existing matrices.
  static int count;
#define FASTMAT_DEBUG_DESTR
#ifdef FASTMAT_DEBUG_DESTR
  int ident;
  static int lastid;
  /// flags if matrix was already defined
  int defined;
  /// flags if matrix was already defined
  int is_defined(void) const {return defined;};
#endif
};

/** @name FastMat operations */
//@{

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Product of matrices. A=B*C
    @author M. Storti
    @param A (output) result of product
    @param B (input) first matrix
    @param C (input) second matrix
*/ 
int FMp(FastMat & A,const FastMat & B,const FastMat & C);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Addition of matrices. A = B + C
    @author M. Storti
    @param A (output) result of sum
    @param B (input) first matrix
    @param C (input) second matrix
*/ 
int FMa(FastMat & A,const FastMat & B,const FastMat & C);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** AXPY operation for fastmat matrices y = a x + y
    @author M. Storti
    @param Y (input/output) matrix to be modified
    @param a (input) scalar
    @param (input) matrix to be added
*/ 
int FMaxpy(FastMat & Y,double a,FastMat & X);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** AXPY operation for fastmat matrices y = a x' + y
    @author M. Storti
    @param Y (input/output) matrix to be modified
    @param a (input) scalar
    @param (input) matrix to be added (transposed)
*/ 
int FMaxpy_t(FastMat & Y,double a,const FastMat & X);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Inverse of a matrix
    @author M. Storti
    @param A (input) the matrix to be inverted
    @param invA (output) the inverted matrix
*/ 
int FMinv(FastMat & invA,FastMat & A);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Transpose of a matrix
    @author M. Storti
    @param A (input) the matrix to be transposed
    @param At (output) the transposed matrix
*/ 
int FMtr(FastMat & At,const FastMat & A);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Determinant of a Fastmat matrix 
    @author M. Storti
    @param A (input) the matrix
    @param det (ouput) the determinant
*/ 
int FMdet(double & det,const FastMat & A);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Converts Newmat matrices to Fastmat matrices
    @author M. Storti
    @param B (input)  matrix in Newmat form
    @param A (output)  matrix in Fastmat form
*/ 
int NM2FM(FastMat & A,const Matrix & B);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Clon of Matlab's kron.
    For given matrices A (nxm) and B (pxq) returns a matrix C (np x
    mq) formed with blocks [A(1,1)*B A(1,2)*B .... ; A(2,1)*B A(2,2)*B
    ; .... A(n,m)*B]
    @author M. Storti
    @param A first matrix argument
    @param B second matrix argument
    @return C=cr/n(A,B)
*/ 
int kron(FastMat & C, FastMat const &A, FastMat const & B);

//@}
#endif

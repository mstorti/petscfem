// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse.h,v 1.11 2001/09/22 14:02:04 mstorti Exp $
#ifndef SEQMAT_H
#define SEQMAT_H

#include <cstdio>

#include <map>
#include <vector>
#include <algorithm>

#include <randomg.h>

using namespace Random;

namespace Sparse {

  class Indx : public vector<int> {
  public:
    Indx(int n,int m,int k=1);
  };
  typedef Indx::iterator IndxIt;
  typedef Indx::const_iterator IndxCIt;

  typedef map<int,double>::iterator VecIt;
  typedef map<int,double>::const_iterator VecCIt;
  typedef pair<int,double> VecP;

  class Mat;
      
  /// A simple sparse vector class (0 indexed). 
  class Vec : public map<int,double> {
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

  };
  
  typedef map<int,Vec>::iterator RowIt;
  typedef map<int,Vec>::const_iterator RowCIt;
  typedef pair<int,Vec> RowP;
  typedef pair<const int,Vec> RowCP;

    // Simple sparse matrix class. 
  class Mat : public map< int, Vec >  {

    /// Dimensions
    int nrows,ncols;
    /// Flag indicating where you can add values past the specified dimensions 
    int grow_m; 
  public:
    friend class Vec;

    /// Constructor from the length
    Mat(int m=0,int n=0) : grow_m(1), nrows(m), ncols(n) {};
    /// Return row dimension
    int rows() {return nrows;};
    /// Return column dimension
    int cols() {return ncols;};
    /// Constructor from another vector
    Mat(const Mat &a) {*this = a;};

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

    /// print elements (sparse version)
    void print(const char *s = NULL);
    /// print elements (full version)
    void print_f(const char *s = NULL);

    /// Resize vectors, truncates elements if greater than this value
    Mat & resize(int m,int n);
    /// Clears all elements
    Mat & clear() {map<int,Vec>::clear(); return *this;}
    /// Set mode if can grow automatically or not
    Mat & grow(int g) { grow_m=g; return *this;};
    /// Flags if the vector is empty or not
    int empty() const;
    /// Number of non null elements
    int size() const;
  };

}
#endif

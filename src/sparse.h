// -*- mode: C++ -*- 
/*__INSERT_LICENSE__*/
// $Id: sparse.h,v 1.5 2001/09/21 15:37:47 mstorti Exp $
#ifndef SEQMAT_H
#define SEQMAT_H

#include <cstdio>

#include <map>
#include <vector>
#include <algorithm>

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
      
  /// A simple sparse vector class (0 indexed). 
  class Vec : public map<int,double> {
    /// Length of the sparse vector
    int len;
    /// Flag indicating where you can add values past the specified length or not 
    int grow_m; 
  public:
    
    /// Constructor from the length
    Vec(int l=0) : grow_m(1) {len=l;};
    /// Return length of the vector 
    int length() {return len;};
    /// Constructor from another vector
    Vec(const Vec &v) {*this = v;};
    /* Insert contents of vector v at position I.
	Length of vector is that of the maximum index in I. 
    */
    // Vec(const Indx &I,const Vec &v);
    /// Get element at specified position
    double get(int j) const;
    /// Get a subvector of elements at position I
    void get(const Indx &I,Vec &V) const; 

    /// Set element at position j
    Vec & set(int j,double v);
    /// Set elements at subvector at position I. w[I] = v
    Vec & set(const Indx & I,const Vec & v);
    /// Set vector to elements of v at position I. w = v[I]
    Vec & set(const Vec & v,const Indx & I);
    /// Set vector elements at I to elements f of v at position J. w[I] = v[J]
    Vec & set(const Indx & I,const Vec & v,const Indx & K);

    /// Scale elements 
    Vec & scale(double c);

    /// Sets w += a * v
    Vec & axpy(double a,const Vec & v);
    /// Sets w[I] = c * w[I] + a * v
    Vec & axpy(double c,const Indx & I,double a,const Vec & v);

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

//      /// Equal operator
//      Vec operator=(Vec v) { *this = v; return *this;}

  };
  
  typedef map<int,Vec>::iterator RowIt;
  typedef map<int,Vec>::const_iterator RowCIt;
  typedef pair<int,Vec> RowP;
  typedef pair<const int,Vec> RowCP;

#if 0
  class IntTrunc {
  private:
    int bound;
  public:
    IntTrunc(int b) : bound(b) {};
    int operator()(pair<int,Vec> j) { return j.first >= bound;}
  };
#endif

    // Simple sparse matrix class. 
  class Mat : public map< int, Vec >  {

    /// Dimensions
    int nrows,ncols;
    /// Flag indicating where you can add values past the specified dimensions 
    int grow_m; 
  public:
    /// Constructor from the length
    Mat(int m=0,int n=0) : grow_m(1), nrows(m), ncols(n) {};
    /// Return row dimension
    int rows() {return nrows;};
    /// Return column dimension
    int cols() {return ncols;};
    /// Constructor from another vector
    Mat(const Mat &a) {*this = a;};

    /// Get element at specified position
    double get(int j,int k) const;
    /// Set element at position j
    Mat & set(int j,int k,double v);

    /// print elements
    void print(const char *s = NULL);

    /// Resize vectors, truncates elements if greater than this value
    Mat & resize(int m,int n);
  };

}
#endif

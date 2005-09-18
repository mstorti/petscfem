// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: utils.h,v 1.19 2005/09/18 20:34:09 mstorti Exp $

#ifndef UTILS_H
#define UTILS_H

#include <set>

#include <petscvec.h>

#include "fem.h"
#include "fstack.h"
#include "util2.h"

/**@name Utils */
//@{

/** Computes the determinant of a newmat matrix. 
    Currently only for n x n with n<=3 matrices. 
    @author M. Storti
    @param A matrix to compute the determinant
    @return deterinant of A
*/ 
double mydet(Matrix A);

/** Computes the absolute value of the normal vector (differential of
    surface).  When integrating over a surface (or a line in 2D) we
    need the differential of surface, which plays the role of the
    differential of volume in 3D which is related with the determinant
    of the Jacobian from master to spatial coordinates. In this case
    we have #ndim# spatial dimensions and #ndim-1# coordinates in the
    master element. The jacobian has dimensions #ndim x (ndim-1)# and
    computing the minors of this matrix we obtain a vector of
    dimension ndim. The norm of this vector is the differential of
    surface. If you need the unit vector normal to the surface, then
    simply normalize this vector. 
    @author M. Storti 
    @param A matrix to compute the determinant 
    @param S vector normal to the surface (with absolute value equal
    to the relative area)
    @return norm of vector formed with the minors of A (norm of normal
    vector in 3D)
*/
double mydetsur(Matrix &A, ColumnVector & S);

/// Overloaded for Newmat matrices.
double mydetsur(FastMat2 &A, FastMat2 &S);

/** Clon of Matlab's kron.
    For given matrices A (nxm) and B (pxq) returns a matrix C (np x
    mq) formed with blocks [A(1,1)*B A(1,2)*B .... ; A(2,1)*B A(2,2)*B
    ; .... A(n,m)*B]
    @author M. Storti
    @param A first matrix argument
    @param B second matrix argument
    @return C=kron(A,B)
*/ 
Matrix kron(const Matrix & A,const Matrix & B);

/** Double random function based on rand().
    Random value uniformly distributed in [0,1].
    @name drand
    @author M. Storti
    @return random value
*/ 
inline double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

/** Integer random number uniformly distributed in [imin,imax]
  @author M. Storti
  @param imin low end  of the range
  @param imin high end of the range
  @return integer random value
*/ 
int irand(int imin,int imax);

/** Integer random number uniformly distributed in [0,n-1]
  @author M. Storti
  @param n number of objects
  @return integer random value
*/ 
int irand(int n);

/** Minimum value of a set of integers.
  @author M. Storti
  @param n number of integer values
  @param ... set of integers
  @return the minimum of these integers
*/ 
int mini(int n,...);

/** Maximum value of a set of integers.
  @author M. Storti
  @param n number of integer values
  @param ... set of integers
  @return the maximum of these integers
*/ 
int maxi(int n,...);

/** Maximum value of a set of doubles.
  @author M. Storti
  @param n number of double values
  @param ... set of doubles
  @return the maximum of these integers
*/ 
double maxd(int n,...);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Reshapes a matrix to be m x n. (Clon of matlab reshape()).
    @author M. Storti
    @param A (input/output) matrix to be reshaped.
    @param m (input) new row dimension 
    @param n (input) new column dimension 
*/ 
int reshape(Matrix &A,int m,int n);

/** Generic algorithm to return a random item from a set.
    @name random\_pop
    @author M. Storti
    @param T class of elements in the set
    @param Tset the set
    @return random element in set
*/
template<typename T> 
T random_pop(set<T> &Tset) {
  typedef typename set<T>::iterator T_iterator;
  T retval;
  int n=Tset.size();
  int j = irand(1,n);
  T_iterator it=Tset.begin();
  for (int k=1; k<j; k++) it++;
  retval = *it;
  Tset.erase(it);
  return retval;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Voids a generic container.
    @author M. Storti
    @param x generic container (set, list, etc...)
*/ 
//#define VOID_IT(x) (x).erase((x).begin(),(x).end())
#define VOID_IT(x) (x).clear()
//@}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int wait_from_console(char *s=NULL);

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Performs the modulo operation, in the sense of number theory. 
    Given integers #k# and #n>0#, we find #m# and #div# such that
    #k = m * div + n# with #0 <= m < n#
    @param k (input) the number to take the modulo
    @param n (input) the modulo
    @param div (input) pointer to an integer where to put the divisor
    @return the result of the modulo operation */ 
inline int modulo(int k, int n, int *div=NULL) {
  int m = k % n;
  int d = k / n;
  if (m<0) {
    m += n;
    d -= 1;
  }
  if (div) *div = d;
  return m;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Performs the modulo operation, in the sense of 
    number theory for doubles. 
    Given doubles #a# and #b>0#, we find #r# and #m# such that
    #a = m * b + r# with #0 <= r < |b|#
    @param k (input) the number to take the modulo
    @param n (input) the modulo
    @param div (input) pointer to an integer where to put the divisor
    @return the result of the modulo operation */ 
inline double modulo(double a, double b,int &m) {
  double r = fmod(a,b);
  // I don't know which is the form to correctly convert
  // here the double `x=(a-r)/b' to an `int'. `x' is almost an 
  // integer, but due to rounding errors the conversion may be
  // to the lowest integer. I then use the 0.5 extra in order to
  // be sure to get the correct number. 
  m = int((a-r)/b+0.5);
  if (r<0) {
    if (b>0) {
      r += b;
      m -= 1;
    } else {
      r -= b;
      m += 1;
    }
  }
  return r;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
inline double modulo(double a, double b) {
  int m;
  return modulo(a,b,m);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Broadcasts a string from the master to the slaves
    @param s (input/output) the string to be broadcasted
    @param master (input) the index of the master
    @param comm (input) the MPI communicator
    @return error code */ 
int string_bcast(string &s,int master,MPI_Comm comm);

#endif

// -*- mode: C++ -*-
/*__INSERT_LICENSE__*/
//$Id: utils.h,v 1.9 2001/08/02 01:54:01 mstorti Exp $
 

#ifndef UTILS_H
#define UTILS_H

#include <set>

#include <vec.h>

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
    need the the differential of surface, which plays the role of the
    differential of volume in 3D which is related with the determinant
    of the Jacobian from master to spatial coordinates. In this case
    we have `ndim' spatial dimensions and `ndim-1' coordinates in the
    master element. The jacobian has dimensions ndim x (ndim-1) and
    computing the minors of this matrix we obtain a vector of
    dimension ndim. The norm of this vector is the differential of
    surface. 
    @author M. Storti 
    @param A matrix to compute the determinant 
    @param S vector normal to the surface
    @return norm of vector formed with the minors of A (norm of normal
    vector in 3D)
*/
double mydetsur(Matrix &A, ColumnVector & S);

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
template<class T> 
T random_pop(set<T> &Tset) {
  T retval;
  int n=Tset.size();
  int j = irand(1,n);
  set<T>::iterator it=Tset.begin();
  for (int k=1; k<j; k++) it++;
  retval = *it;
  Tset.erase(it);
  return retval;
}

/** Voids a generic container.
    @author M. Storti
    @param x generic container (set, list, etc...)
*/ 
//#define VOID_IT(x) (x).erase((x).begin(),(x).end())
#define VOID_IT(x) (x).clear()
//@}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int wait_from_console(char *s=NULL);

#endif

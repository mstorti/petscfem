// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: cloud2.h,v 1.2 2003/02/27 03:32:41 mstorti Exp $
#ifndef PETSCFEM_CLOUD2_H
#define PETSCFEM_CLOUD2_H

#include <src/fastmat2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Given certain number of points #x# in dimension #ndim# it gives
    the coefficients #w# of the approximation to the a certain set of
    derivatices #nderiv# derivative at point #x0# with a least square
    approximation to a polynomial of order #npol#.  Basically, you
    should have #nx>npol>=nderiv#, and the precision of the
    approximation is #O(npol-nderiv)#.  The computations are computed
    in a FastMat2 loop. */
class BasicCloud {
public:
  /** Dtor. */
  virtual ~BasicCloud()=0;
  /** Initialize the chain. Execute first to all calls to #coef()#
      @param ndim (input) Spatial dimension
      @param nx (input) Number of points 
      @param nderiv (input) How much derivatives we will get.
      @param derivs (input) Array of #ndim*nderiv# integers indicating which
      derivatives we want. For instance in 2D, #nderiv=2# and 
      #derivs={1,1,2,3}# means obtaning approximation to #phi_{xy}# and
      #phi_{xxyyy}#. 
      @param npol (input) Array of #ndim# integers indicating the order of the 
      polynomial to be used, for instance in 2D #npol={2,3}# means to take polynomials
      with order up to #x^2*y^3#. */ 
  virtual void init(int ndim, int nx, int nderiv,const int *derivs, const int *npol)=0;
  /** Compute the coefficients for a certain set of coordinates. 
      @param x (input) coordinates of points (size #nx*ndim#)
      @param w (output) coefficients of approximation (size #nx*nderiv#)
      @param x0 (input) center where the derivative is approximated. (size #ndim#) */
  virtual void coef(FastMat2 &x, FastMat2 &w,FastMat2 &x0)=0;
  /** Compute the coefficients for a certain set of coordinates. Same as 
      #coef(x,w,x0)# but assumes #x0=0#. 
      @param x (input) coordinates of points (size #nx*ndim#)
      @param w (output) coefficients of approximation (size #nx*nderiv#) */
  void coef(FastMat2 &x, FastMat2 &w) { /* fixme:= code */ }
  /** Clean-up function. */ 
  virtual void clear()=0;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Cloud2 : public BasicCloud {
private:
  BasicCloud *ptr;
public:
  /** Ctor. */
  Cloud2();
  /** Dtor. */
  ~Cloud2();
  void init(int ndim, int nx, int nderiv,const int *derivs, const int *npol);
  void coef(FastMat2 &x, FastMat2 &w, FastMat2 &x0);
  void clear();
};

#endif

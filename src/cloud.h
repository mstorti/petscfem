// -*- mode: c++ -*-
//__INSERT_LICENSE__
// $Id: cloud.h,v 1.2 2003/02/25 03:14:04 mstorti Exp $
#ifndef PETSCFEM_CLOUD_H
#define PETSCFEM_CLOUD_H

#include <src/fastmat2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
/** Given certain number of points #x# (1D) it gives the coefficients
    #w# of the approximation to the #nderiv# derivative at point #x0#
    with a least square approximation to a polynomial of order #npol#.
    Basically, you should have #nx>npol>=nderiv#, and the precision of
    the approximation is #O(npol-nderiv)#.  The computations are
    computed in a FastMat2 loop. */
class Cloud {
private:
  /// Order of derivative
  int nderiv;
  /// Order of polynomial
  int npol;
  /// Number of points 
  int nx;
  /// internal 
  double nderiv_fact;
  FastMat2 A, xi, H, iH, AA;
public:
  /** Initialize the chain. Execute first to all calls to #coef()#
      @param nx (input) Number of points 
      @param nderiv (input) order of derivative
      @param mpol (input) order of polynomial to be used. */ 
  void init(int nx, int nderiv,int npol);
  /** Compute the coefficients for a certain set of coordinates. 
      @param x (input) coordinates of points
      @param w (output) coefficients of approximation
      @param x0 (input) center where the derivative is approximated. */
  void coef(FastMat2 &x, FastMat2 &w,double x0=0.);
  /** Clean-up function. */ 
  void clear();
};

#endif

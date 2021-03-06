//__INSERT_LICENSE__
//$Id: secant.cpp,v 1.5 2003/07/02 23:22:19 mstorti Exp $

#include <cstdio>
#include <cmath>
#include "secant.h"

// Solves the equation
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "double Secant::sol(void)" 
double Secant::sol(void) {
  assert(epsilon>0);
  // The sequence order is x1 -> x0 -> x2
  double f0,f1,x1,x2,f2; 
  x1 = x0+epsilon;   x2=x0;
  f0 = residual(x0); f2=f0;
  f1 = residual(x1);
  int j;
  for (j=0; j<maxits; j++) {
    double slope = (f1-f0)/(x1-x0);
    x2=x0-f0/slope;
    f2=residual(x2);
    if (std::abs(f2)<tol) {
      break;
    }
    x1=x0;
    f1=f0;
    x0=x2;
    f0=f2;
  }
  its = j;
  f=f2;
  x0=x2;			// Store last solution for
				// initialization in next call
  return x2;
}

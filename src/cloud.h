//__INSERT_LICENSE__
// $Id: cloud.h,v 1.1 2003/02/25 00:06:11 mstorti Exp $
#ifndef PETSCFEM_CLOUD_H
#define PETSCFEM_CLOUD_H

#include <src/fastmat2.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Cloud {
private:
  int nderiv, npol, nx;
  double nderiv_fact;
  FastMat2 A, xi, H, iH, AA;
public:
  void init(int nx_a, int nderiv_a,int npol_a);
  void coef(FastMat2 &x, FastMat2 &w,double x0=0.);
  void clear();
};

#endif

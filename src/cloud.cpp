//__INSERT_LICENSE__
// $Id: cloud.cpp,v 1.1 2003/02/25 00:06:11 mstorti Exp $
#include <cmath>
#include <src/cloud.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static double my_power(double x,void* p) {
  return pow(x,*(double *)p);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud::init(int nx_a, int nderiv_a,int npol_a) {
  nderiv = nderiv_a;
  npol = npol_a;
  nx = nx_a;
  A.resize(2,nx,npol+1);
  AA.resize(2,nx,npol+1);
  xi.resize(1,nx);
  H.resize(2,npol+1,npol+1);
  iH.resize(2,npol+1,npol+1);
  nderiv_fact = 1.;
  for (int k=2; k<=nderiv; k++)
    nderiv_fact *= double(k);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud::coef(FastMat2 &x, FastMat2 &w,double x0=0.) {
  A.ir(2,1).set(1.);
  xi.set(x).add(-x0);
  double h = xi.max_abs_all();
  xi.scale(1./h);
  for (int k=1; k<=npol; k++) {
    double exp = double(k);
    A.ir(2,k+1).set(xi).fun(&my_power,&exp);
  }
  A.rs();
  H.prod(A,A,-1,1,-1,2);
  iH.inv(H);
  AA.prod(A,iH,1,-1,2,-1);
  AA.ir(2,nderiv+1);
  w.set(AA).scale(nderiv_fact/pow(h,nderiv));
  AA.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud::clear() { 
  A.clear(); 
  xi.clear(); 
  H.clear();
  iH.clear();
}

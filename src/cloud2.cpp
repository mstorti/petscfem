//__INSERT_LICENSE__
// $Id: cloud2.cpp,v 1.4 2003/02/28 23:51:05 mstorti Exp $
#include <cmath>
#include <src/util2.h>
#include <src/dvector.h>
#include <src/dvector2.h>
#include <src/cloud2.h>
#include <src/fm2temp.h>

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class Cloud2FastMat2 : public BasicCloud {
private:
  /// Space dimension
  int ndim;
  /// Number of points 
  int nx;
  /// Number of derivative coefficients to be computed
  int nderiv;
  /** Array of #nderiv*ndim# integers indicating which
      derivatives we want. */
  dvector<int> derivs;
  /** Array of #nderiv*ndim# integers indicating which
      derivatives we want. */
  dvector<int> coef0;
  /// Order of polynomial
  dvector<int> npol;
  /** Indices of the column corresponding to the
      requested derivatives */
  dvector<int> deriv_indx;
  /// Exponent in scaling factor #h^sum_j(k_j)#
  dvector<int> h_fac_expo;
  /// Factorial factor #prod_j k_j!#
  dvector<double> n_fact;
  /// Number of terms in polynomial
  int n_pol_term;
  /// To store temporaries
  FastMat2Tmp tmp;
  /// Aux matrices
  FastMat2 A,AA,xi,H,iH,x00;
public:
  void init(int ndim, int nx, int nderiv,const int *derivs, const int *npol);
  void coef(FastMat2 &x, FastMat2 &w,FastMat2 &x0);
  void coef(FastMat2 &x, FastMat2 &w);
  double cond();
  void clear();
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
BasicCloud::~BasicCloud() {}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Cloud2::Cloud2() : ptr(NULL) { ptr = new Cloud2FastMat2; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
Cloud2::~Cloud2() { delete ptr; }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2::init(int ndim, int nx, int nderiv,const int *derivs, const int *npol) {
  ptr->init(ndim,nx,nderiv,derivs,npol);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2::coef(FastMat2 &x, FastMat2 &w,FastMat2 &x0) { ptr->coef(x,w,x0); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2::coef(FastMat2 &x, FastMat2 &w) { ptr->coef(x,w); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Cloud2::cond() { return ptr->cond(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2::clear() { ptr->clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
static double my_power(double x,void* p) {
  int e = *(int *)p;
  return int_pow(x,e);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
int facto(int n) {
  assert(n>=0);
  int fac = 1;
  for (int j=1; j<=n; j++) fac *= j;
  return fac;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2FastMat2::init(int ndim_a, int nx_a, 
			  int nderiv_a,const int *derivs_a, const int *npol_a) {
  ndim = ndim_a;
  nx = nx_a;
  nderiv = nderiv_a;

  npol.a_resize(1,ndim);
  npol.set(npol_a);

  derivs.a_resize(2,nderiv,ndim);
  derivs.set(derivs_a);
  
  deriv_indx.resize(nderiv);
  h_fac_expo.resize(nderiv);
  n_fact.resize(nderiv);
  for (int j=0; j<nderiv; j++) {
    int index = derivs.e(j,ndim-1);
    for (int l=ndim-2; l>=0; l--) {
      assert(derivs.e(j,l)<=npol.e(l));
      index = index * (npol.e(l)+1) + derivs.e(j,l);
    }

    int f = 1;
    int ee = 0;
    for (int l=0; l<ndim; l++) {
      f *= facto(derivs.e(j,l));
      ee += derivs.e(j,l);
    }
    deriv_indx.e(j) = index;
    n_fact.e(j) = double(f);
    h_fac_expo.e(j) = ee;
  }

  // Total number of terms is the product of all the (orders+1)
  n_pol_term=1;
  for (int j=0; j<ndim; j++) n_pol_term *= npol.e(j)+1;

  A.resize(2,nx,n_pol_term);
  AA.resize(2,nx,n_pol_term);
  xi.resize(2,nx,ndim);
  H.resize(2,n_pol_term,n_pol_term);
  iH.resize(2,n_pol_term,n_pol_term);
  tmp(0).resize(1,nx).set(1.);
  tmp(3).resize(1,nx);
  x00.resize(1,ndim).set(0.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2FastMat2::coef(FastMat2 &x, FastMat2 &w) {
  coef(x,w,x00);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2FastMat2::coef(FastMat2 &x, FastMat2 &w,FastMat2 &x0) {
  // xi = (x-x0)/h
  tmp(1).prod(tmp(0),x0,1,2);
  xi.set(x).rest(tmp(1));
  double h = xi.max_abs_all();
  xi.scale(1./h);

  dvector<int> expo;
  expo.a_resize(1,ndim).set(int(0));
  int ideriv = 0;
  while (1) {
    tmp(2).set(xi);
    tmp(3).set(1.);
    for (int j=0; j<ndim; j++) {
      tmp(2).ir(2,j+1).fun(&my_power,&expo.e(j));
      tmp(3).mult(tmp(2));
    }
    tmp(2).rs();
    A.ir(2,ideriv+1).set(tmp(3));
    int j;
    for (j=0; j<ndim; j++) { // Increment
      expo.e(j)++;
      if (expo.e(j)<=npol.e(j)) break;
      expo.e(j)=0;
    }
    ideriv++;
    if (j==ndim) break;
  }
  A.rs();
  H.prod(A,A,-1,1,-1,2);
  iH.inv(H);
  AA.prod(A,iH,1,-1,2,-1);
  for (int j=0; j<nderiv; j++) {
    AA.ir(2,deriv_indx.e(j)+1);
    w.ir(2,j+1).set(AA).scale(n_fact.e(j)/int_pow(h,h_fac_expo.e(j)));
  }
  AA.rs();
  w.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double Cloud2FastMat2::cond() { 
  return  H.norm_p_all(1.)*iH.norm_p_all(1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void Cloud2FastMat2::clear() { 
  A.clear(); 
  xi.clear(); 
  H.clear();
  iH.clear();
  derivs.clear();
  npol.clear();
  tmp.clear();
  x00.clear();
}

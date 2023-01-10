//__INSERT_LICENSE__
//$Id: linhff.cpp,v 1.9 2003/11/25 01:13:36 mstorti Exp $
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "advective.h"
#include "genload.h"

// Regularized version of the abs function
static double regabs(double x,double delta=1e-4) {
  return  (fabs(x)<1e-6? 1.0 : x/tanh(x));
}

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double regmin(double a,double b,double delta=1e-4) {
  return 0.5*(a+b)-0.5*regabs(a-b,delta);
}
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double regmax(double a,double b,double delta=1e-4) {
  return 0.5*(a+b)+0.5*regabs(a-b,delta);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double fluxfun(double DV) {
  double R0=100,Rinf=0.01,DV0=1,delta=0.1;
  return regmax(DV/R0,(DV-DV0)/Rinf,delta);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {
  static int k=0;
  k++;
#if 0
  dU.set(uout).minus(uin);
  h->prod(flux,dU);
  h->jac(jacin);
  jacout.set(jacin);
  jacout.scale(-1.);
  if (k>10) {
    FMSHV(uin);
    FMSHV(uout);
    FMSHV(flux);
    FMSHV(jacin);
    FMSHV(jacout);
    exit(0);
  }
  s->add(flux);
#else
  // printf("in linhff\n");
  double
    *uinp = uin.data(),
    *uoutp = uout.data(),
    *fluxp = flux.data(),
    *jacinp = jacin.data(),
    *jacoutp = jacout.data();
  // double hfilm=2.0;
  double DV = (*uoutp-*uinp);
  double epsln = 1e-5;
  *fluxp = fluxfun(DV);
  double hfilm =(fluxfun(DV+epsln)-fluxfun(DV-epsln))/(2*epsln);
  *jacinp = hfilm;
  *jacoutp = -hfilm;
  // if (k>10) {  }
  // exit(0);
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin) {
  dU.set(uin).scale(-1.);
  h->prod(flux,dU);
  h->jac(jacin);
  s->add(flux);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::HFull::element_hook(ElementIterator &element) {
  const double * hf = l->elemset->prop_array(element,l->hfilm_coeff_prop);
  HH.set(hf);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::HFull::init() {
  HH.resize(2,l->ndof,l->ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::HFull::jac(FastMat2 &A) {
  A.set(HH);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::SFull::element_hook(ElementIterator &element) {
  const double * s = l->elemset->prop_array(element,l->hfilm_source_prop);
  SS.set(s);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::init() {
  elemset->elem_params(nel,ndof,nelprops);
  // Read hfilm coefficients. 
  //o _T: double[var_len]
  //  _N: hfilm_coeff _D: no default  _DOC: 
  // Defines coeffcients for the film flux function. May be 
  //  #var_len=0#  (no $\Delta T$ driven load) or
  //  #var_len=ndof*ndof#  a full matrix of relating the flux with
  // $\Delta !U$. 
  //  _END
  elemset->get_prop(hfilm_coeff_prop,"hfilm_coeff");
  if (hfilm_coeff_prop.length == ndof*ndof) {
    h= new HFull(this);
  } else if (hfilm_coeff_prop.length == 0) {
    h= new HNull(this);
  } else {
    PETSCFEM_ERROR("Not valid size of hfilm_coeff: %d, ndof: %d\n",
		   hfilm_coeff_prop.length,ndof);
  }

  // Read source term for generic load elemset. 
  //o _T: double[var_len]
  //  _N: hfilm_source _D: no default  _DOC: 
  // Defines constant source term for the generic load on
  // surfaces. May be of length 0 (null load) or  #ndof# 
  // which represents a geven load per field. 
  //  _END
  elemset->get_prop(hfilm_source_prop,"hfilm_source");
  if (hfilm_source_prop.length == ndof) {
    s= new SFull(this);
  } else if (hfilm_source_prop.length == 0) {
    s= new SNull(this);
  } else {
    PETSCFEM_ERROR("Not valid size of hfilm_source: %d, ndof: %d\n",
		   hfilm_source_prop.length,ndof);
  }

  dU.resize(1,ndof);
  h->init();
  s->init();
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
LinearHFilmFun::~LinearHFilmFun() {
  delete h;
  delete s;
}

//__INSERT_LICENSE__
//$Id: linhff.cpp,v 1.6 2001/12/20 21:58:51 mstorti Exp $
 
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "advective.h"
#include "genload.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {

  dU.set(uout).rest(uin);
  h->prod(flux,dU);
  h->jac(jacin);
  jacout.set(jacin);
  jacout.scale(-1.);
  s->add(flux);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &flux,FastMat2 &jacin) {
  jacin.set(0.);
  flux.set(0.);
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
#undef __FUNC__
#define __FUNC__ "void LinearHFilmFun::init()" 
void LinearHFilmFun::init() {
  elemset->elem_params(nel,ndof,nelprops);
  // Read hfilm coefficients. 
  //o _T: double[var_len]
  //  _N: hfilm_coeff _D: no default  _DOC: 
  // Defines coeffcients for the film flux function. May be 
  // \verb+var_len=0+ (no $\Delta T$ driven load) or
  // \verb+var_len=ndof*ndof+ a full matrix of relating the flux with
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
  // surfaces. May be of length 0 (null load) or \verb+ndof+
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

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void LinearHFilmFun::q(...)"
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "void LinearHFilmFun::q(...)"
void LinearHFilmFun::init(void) {
}
#endif

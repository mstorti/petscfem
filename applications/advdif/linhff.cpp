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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {
  // FLAG is for doing the initialization just once
  // USE_ELYZER_FILM is to flag if the special nonlinear functions
  // must be taken
#if 0  
  if (!fluxfun.flag) {
    fluxfun.flag=1;
    int ierr;
    int &uef = fluxfun.use_elyzer_film;
    uef = get_entry_d(elemset,"use_elyzer_film");
    // TGETOPTDEF(elemset->thash,int,use_elyzer_film,0);
    // fluxfun.use_elyzer_film = use_elyzer_film;
    if (uef) fluxfun.init(elemset);
    if (!MY_RANK) printf("elemset %p use_elyzer_film %d\n",elemset,uef);
  }
#endif
  if (fluxfun_table[elemset]) {
    // A fluxfun may have been set for this elemset in a hook.
    // If this is so move it here, to the elemset
    TRACE("Moving fluxfun ptr");
    fluxfunp = move(fluxfun_table[elemset]);
  }
  
  if (!fluxfunp.get()) {
    // Use the normal linear functions
    dU.set(uout).minus(uin);
    h->prod(flux,dU);
    h->jac(jacin);
    jacout.set(jacin);
    jacout.scale(-1.);
#if 0    
    if (0) {
      FMSHV(uin);
      FMSHV(uout);
      FMSHV(flux);
      FMSHV(jacin);
      FMSHV(jacout);
      exit(0);
    }
#endif
    s->add(flux);
  } else {
    static int flag=0;
    if (!flag) {
      flag = 1;
      TRACE("Using fluxfun ptr");
    }
    // Use the functions provided by the user
    // We dont use FastMat2 so we get the pointers to the
    // internal data
    double
      *uinp = uin.data(),
      *uoutp = uout.data(),
      *fluxp = flux.data(),
      *jacinp = jacin.data(),
      *jacoutp = jacout.data();
    // Difference potential about this film
    double DV = (*uoutp-*uinp);
    // Small increment to take the Jacobian by finite differences
    double epsln = 1e-5;

    // Call the fluxfun function to get the flux
    *fluxp = fluxfunp->fun(DV,1);
    // Compute the Jacobian by finite differences
    double hfilm =(fluxfunp->fun(DV+epsln)-fluxfunp->fun(DV-epsln))/(2*epsln);
    // Set the Jacobians
    *jacinp = hfilm;
    *jacoutp = -hfilm;
  }
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

map<NewElemset*,unique_ptr<fluxfun_t>> LinearHFilmFun::fluxfun_table;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::init() {
  static int call=0;
  printf("In LinearHFilmFun::init, call %d\n",call++);
  elemset->elem_params(nel,ndof,nelprops);
  // Read hfilm coefficients. 
  //o _T: double[var_len]
  //  _N: hfilm_coeff _D: no default  _DOC: 
  // Defines coeffcients for the film flux function. May be 
  //  #var_len=0#  (no $\Delta T$ driven load) or
  //  #var_len=ndof*ndof#  a full matrix of relating the flux with
  // $\Delta !U$. 
  //  _END
  if (!h.get()) {
    elemset->get_prop(hfilm_coeff_prop,"hfilm_coeff");
    if (hfilm_coeff_prop.length == ndof*ndof) {
      h.reset(new HFull(this));
    } else if (hfilm_coeff_prop.length == 0) {
      h.reset(new HNull(this));
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
      s.reset(new SFull(this));
    } else if (hfilm_source_prop.length == 0) {
      s.reset(new SNull(this));
    } else {
      PETSCFEM_ERROR("Not valid size of hfilm_source: %d, ndof: %d\n",
                     hfilm_source_prop.length,ndof);
    }
  }
  
  dU.resize(1,ndof);
  h->init();
  s->init();
  // Just set the entry in the table
  fluxfun_table[elemset];
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
LinearHFilmFun::~LinearHFilmFun() {
  // delete h;
  // delete s;
}

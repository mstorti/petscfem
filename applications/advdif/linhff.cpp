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

static int VRBS=0;

// Regularized version of the abs function
static double regabs(double x,double delta=1e-4) {
  double y = (fabs(x)<1e-6? 1.0 : x/tanh(x));
  if (0 && VRBS) printf("x %g, y %g, delta %g\n",x,y,delta);
  return y;
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
double fluxfun_t::fun(double DV) {
  double
    aDV=fabs(DV),
    sig=(DV>0? 1 : -1),
    DV0 = (sig>0? DV0p : DV0m),
    flux = regmax(aDV/R0,(aDV-DV0)/Rinf,delta);
  if (0 && VRBS) printf("aDV %g, sig %g, DV0 %g, flux %g\n",aDV,sig,DV0,flux);
  return sig*flux;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double get_entry_d(NewElemset *e,const char *name) {
  const char *s=NULL;
  e->get_entry(name,s);
  // printf("blabla %s\n",s);
  PETSCFEM_ASSERT(s!=NULL,"not found entry %s!!",name);
  return stod(s);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void fluxfun_t::init(NewElemset *e) {
  int ierr;
  printf("name %s\n",e->name());
  // const char *s;
  // e->get_entry("blabla",s);
  // printf("blabla %s\n",s);
  // PETSCFEM_ASSERT0(s!=NULL,"not found entry!!");
  // double blabla = stod(s);
#define GET_ENTRY_D(name) name = get_entry_d(e,#name)
  GET_ENTRY_D(R0);
  GET_ENTRY_D(Rinf);
  GET_ENTRY_D(DV0p);
  GET_ENTRY_D(DV0m);
  GET_ENTRY_D(delta);
  // TGETOPTDEF_ND(GLOBAL_OPTIONS,double,R0,NAN);
  // TGETOPTDEF_ND(GLOBAL_OPTIONS,double,Rinf,NAN);
  // TGETOPTDEF_ND(GLOBAL_OPTIONS,double,DV0,NAN);
  // TGETOPTDEF_ND(GLOBAL_OPTIONS,double,delta,NAN);
  if (!MY_RANK) 
    printf("USER FLUXFUN initialized: R0 %g, Rinf %g, "
           "DV0p %g, DV0m %g, delta %g\n",
           R0,Rinf,DV0p,DV0m,delta);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {
  // FLAG is for doing the initialization just once
  // USE_ELYZER_FILM is to flag if the special nonlinear functions
  // must be taken
  static int flag=0,use_elyzer_film=0;
  if (!fluxfun.flag) {
    fluxfun.flag=1;
    int ierr;
    TGETOPTDEF_ND(GLOBAL_OPTIONS,int,use_elyzer_film,0);
    if (use_elyzer_film) fluxfun.init(elemset);
  }

  if (use_elyzer_film==0) {
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

    static int cnt=0; cnt++;
    if (cnt>2000) { 
      int N=1000;
      double a=-1,b=1;
      VRBS = 1;
      for (int j=0; j<N; j++) {
        double
          x = a+double(j)/N*(b-a),
          y = fluxfun.fun(DV);
        printf("%g %g\n",x,y);
      }
      VRBS = 0;
      exit(0);
    }
    
    // Call the function to get the flux
    *fluxp = fluxfun.fun(DV);
    // Compute the Jacobian by finite differences
    double hfilm =(fluxfun.fun(DV+epsln)-fluxfun.fun(DV-epsln))/(2*epsln);
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

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

#if 0
// These functions are commented out because otheriwse the
// compiler complains about non used function

// Regularized version of the abs function
static double regabs(double x,double delta=1e-4) {
  double ax = x/delta,
    y = delta*(fabs(ax)<1e-6? 1.0 : ax/tanh(ax));
  if (0 && VRBS) printf("x %g, y %g, delta %g\n",x,y,delta);
  return y;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double regmin(double a,double b,double delta=1e-4) {
  return 0.5*(a+b)-0.5*regabs(a-b,delta);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double regmax(double a,double b,double delta=1e-4) {
  return 0.5*(a+b)+0.5*regabs(a-b,delta);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
static double regpos2(double x,double delta) {
  double xx = x-delta;
  return 0.5*(xx+pf_regabs(xx,delta));
}
#endif

// This global variable allows to set the Rinf from a hook
double FLUXFUN_H2_RINF=NAN;
lhff_info_t LHH_INFO;

#if 0
//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double fluxfun_h2_t::fun(double DV) {
  // Store the last value so that report when the value changes
  static double Rinf_last = NAN;
  // If the user has set the global value in a hook copy on the used value
  if (!ISNAN(FLUXFUN_H2_RINF)) Rinf = FLUXFUN_H2_RINF;
  // Report if Rinf value was changed
  if (!MY_RANK && Rinf!=Rinf_last) {
    printf("Changed Rinf %f -> %f\n",Rinf_last,Rinf);
    Rinf_last = Rinf;
  }
#if 0
  // OLD VERSION
  double
    aDV=fabs(DV),
    sig=(DV>0? 1 : -1),
    DV0 = (sig>0? DV0p : DV0m),
    flux = regmax(aDV/R0,(aDV-DV0)/Rinf,delta);
  flux *= sig;
  if (0 && VRBS && rand()%1000==0)
    printf("aDV %g, sig %g, DV0 %g, flux %g\n",aDV,sig,DV0,flux);
  return flux;
#else
  // const double Z=0.00819266;
  // const double Z=8.2557e-03;
  // FOR RINF=1e-4
  // Z=0        => G=0.248333
  // Z=0.248333 => G=0.134514
  // Z=0.382847 => G=0.0728718
  // Z=0.455718 => G=0.0394622
  // NEW VERSION TO BE USED WITH ITERATIVE PENALIZATION
  const double delta=0.01;
  auto &epg = LHH_INFO.ELEMPG;
  auto &I = LHH_INFO.table[epg];
  I.gfun = DV-DV0p;
  I.gfunm = -DV-DV0m;
  // printf("in linhff: elem %d ipg %d, ZCURRENT %g, GFUN %g, ZCURRENTM %g, GFUNM %g\n",
  //        epg.first,epg.second,I.zcurrent,I.gfun,I.zcurrentm,I.gfunm);
  double fluxp = regpos2(I.zcurrent+I.gfun,delta)/Rinf;
  double fluxm = -regpos2(I.zcurrentm+I.gfunm,delta)/Rinf;
  return fluxp+fluxm;
#endif
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
void fluxfun_h2_t::init(NewElemset *e) {
  int ierr;
  // printf("name %s\n",e->name());
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
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void LinearHFilmFun::q(FastMat2 &uin,FastMat2 &uout,FastMat2 &flux,
		       FastMat2 &jacin,FastMat2 &jacout) {
#if 0
  // FLAG is for doing the initialization just once
  // USE_ELYZER_FILM is to flag if the special nonlinear functions
  // must be taken
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

#if 0
    auto &f = fluxfun;
    int N=1000;
    double a=0,b=1;
    double delta=0.01;
    for (int j=0; j<N; j++) {
      double
        x = a+double(j)/N*(b-a),
        yflux = f.fun(x),
        ynew = regpos2(x-f.DV0p,delta)/f.Rinf;
      printf("x %g yflux %g ynew %g\n",x,yflux,ynew);
    }
    exit(0);
#endif
    
    // Call the function to get the flux
    *fluxp = fluxfunp->fun(DV);
    auto &epg = LHH_INFO.ELEMPG;
    auto &I = LHH_INFO.table[epg];
    I.glast = I.gfun;
    I.glastm = I.gfunm;
    I.flux = *fluxp;
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
  // Just set the entry in the table
  fluxfun_table[e];
}  

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
LinearHFilmFun::~LinearHFilmFun() {
  delete h;
  delete s;
}

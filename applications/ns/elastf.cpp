//__INSERT_LICENSE__
//$Id: elastf.cpp,v 1.1.2.1 2001/10/29 14:34:41 mstorti Exp $

#include <fem.h>
#include <utils.h>
#include <readmesh.h>
#include <getprop.h>
#include <fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

#define ICONE(j,k) (icone[nel*(j)+(k)]) 

extern "C" {
  void initf_(int *,int *,int *,int *,int *,int *,int *,
	       int *,int *,int *,int *,double *,
	       double *,double *,double *,
	       double *,int *,int *,int *);

  void frin30_(double *cartd,double *dvolu,double *ehist,double *elcod,
	       double *eldis,double *emass,double *epmtx,double *gpcod,
	       int *lnods,double *props,double *rmat1,double *shape,
	       double *stran,double *strsg,double *stra0,double *strs0,
	       double *thick,double *bmsig,double *bmatx,double *desig,
	       double *dmatx,double *dstra,double *presg,double *sgtot,
	       double *sigma,double *tstra,double *xjacm,double *tenod,
	       double *dteno,double *pwoel,double *preal,double *tgapl,
	       double *tenoi,double *bsbar,double *xjacn,double *fpchl);

  void stif30_(double *cartd,double *dvolu,double *ehist,
	       double *elcod,double *eldis,double *epmtx,
	       double *gpcod,double *props,double *shape,
	       double *strsg,double *estif,double *hstif,
	       double *emass,double *bmatx,double *dmatx,
	       double *sigma,double *xjacm,double *tenod,
	       double *bsbar,double *dteno,double *xjacn);

}

void elasticity_f::init() {

  int ierr;

  //o Young modulus
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);

  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,yield_strength,0.);
  assert(yield_strength>=0.);

  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;
  // Dos opciones para imprimir
  // printf("rec_Dt: %d\n",rec_Dt);
  // SHV(rec_Dt);
  assert(!(rec_Dt>0. && rho==0.));
  assert(ndof==ndim);

  nen = nel*ndof;
  dshapex.resize(3,ndim,nel,npg);
  dshapex_f.resize(3,npg,nel,ndim); // aux for index inversion

  wpgdet_v.resize(1,npg);
  Jaco.resize(2,ndim,ndim);
  iJaco.resize(2,ndim,ndim);
  int nhist = 36, nstr1 = 4, nkost = 16, nprop = 100,
    nstre = 3, nstrs = 4, nfpch = 1;

  ehist_f.resize(2,npg,nhist).set(0.);
  ehist_m.resize(2,nhist,npg).set(0.);
  emass_m.resize(1,nen*nen).set(0.);
  epmtx_m.resize(1,nkost*npg).set(0.);
  gpcod_m.resize(2,ndim,npg).set(0.);
  gpcod_f.resize(2,npg,ndim);
  props_m.resize(1,nprop).set(0.);
  rmat1_m.resize(2,ndim,ndim).set(0.);
  desig_m.resize(1,nstr1).set(0.);
  dstra_m.resize(1,nstr1).set(0.);
  presg_m.resize(1,nstr1).set(0.);
  sgtot_m.resize(1,nstr1).set(0.);
  sigma_m.resize(1,nstr1).set(0.);
  tstra_m.resize(1,nstr1).set(0.);
  tenod_m.resize(1,nel).set(0.);
  dteno_m.resize(1,nel).set(0.);
  pwoel_m.resize(1,nel).set(0.);
  tgapl_m.resize(1,nel).set(0.);
  preal_m.resize(1,nel).set(0.);
  tenoi_m.resize(1,nel).set(0.);

  bsbar_m.resize(3,nstr1,nen,npg).set(0.);
  bsbar_f.resize(3,npg,nen,nstr1).set(0.);

#define MAKE_DIM(name,d1,d2)         \
  name##_m.resize(2,d1,d2).set(0.);  \
  name##_f.resize(2,d2,d1).set(0.)
  
  MAKE_DIM(stran,nstr1,npg);
  MAKE_DIM(strsg,nstr1,npg);
  MAKE_DIM(stra0,nstr1,npg);
  MAKE_DIM(strs0,nstr1,npg);
  MAKE_DIM(bmatx,nstre,nen);
  MAKE_DIM(dmatx,nstrs,nstrs);
  MAKE_DIM(xjacm,ndim,ndim);
  MAKE_DIM(xjacn,ndim,ndim);
  MAKE_DIM(fpchl,nfpch,nel);
  MAKE_DIM(hstif,nen,nel);

  props_m.setel(0.,   1); // Isotropic
  props_m.setel(20.,  2); // Plastic
  props_m.setel(1.,   3); // No temperature dependent
  props_m.setel(0.,   4); // Standard
  props_m.setel(3.,   59); // Solidification
  props_m.setel(32.,  36); // von Mises
  props_m.setel(1.,   35); // von Mises version
  props_m.setel(32.,  52); // von Mises
  props_m.setel(1.,   51); // von Mises version
  
  double *props = props_m.storage_begin();
  int kdyna = !glob_param->steady;

  initf_(&ndim,&nel,&npg,&nen,
	  &nhist,&nstr1,&nkost,&nprop,
	  &nstre,&nstrs,&nfpch,props,&E,&nu,
	  &yield_strength,
	  &glob_param->Dt,&kdyna,&glob_param->inwt,
	  &ndof);
  
  thickness = 1.;

}

void elasticity_f::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  
  FastMat2 &shape_m = shape; // Change name to `shape'
  double *cartd,*dvolu,*ehist,*elcod,*eldis,*emass,*epmtx,
    *gpcod,*props,*rmat1,*shape,*stran,*strsg,
    *stra0,*strs0,*thick,
    *bmsig,*bmatx,*desig,*dmatx,*dstra,*presg,*sgtot,
    *sigma,*tstra,*xjacm,*tenod,*dteno,
    *pwoel,*preal,*tgapl,*tenoi,*bsbar,*xjacn,*fpchl,
    *estif,*hstif;
  int *lnods;

  for (int ipg=0; ipg<npg; ipg++) {
    
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco = Jaco.det();
    if (detJaco <= 0.) {
      printf("Jacobian of element %d is negative or null\n"
	     " Jacobian: %f\n",elem,detJaco);
      PetscFinalize();
      exit(0);
    }
    wpgdet_v.setel(detJaco*wpg.get(ipg+1),ipg+1);
    iJaco.inv(Jaco);
    dshapex.ir(3,ipg+1).prod(iJaco,dshapexi,1,-1,-1,2);
  } // end of loop over Gauss points 

  gpcod_m.prod(shape_m,xloc,-1,2,-1,1);
  dshapex.rs();

  // Indexes in C and Fortran are inverted
  dshapex_f.ctr(dshapex,3,2,1);
  ehist_f.ctr(ehist_m,2,1);
  gpcod_f.ctr(gpcod_m,2,1);
  shape_f.ctr(shape_m,2,1);
  stran_f.ctr(stran_m,2,1);
  strsg_f.ctr(strsg_m,2,1);
  stra0_f.ctr(stra0_m,2,1);
  strs0_f.ctr(strs0_m,2,1);
  bmatx_f.ctr(bmatx_m,2,1);
  dmatx_f.ctr(dmatx_m,2,1);
  xjacm_f.ctr(xjacm_m,2,1);
  xjacn_f.ctr(xjacn_m,2,1);
  bsbar_f.ctr(bsbar_m,3,2,1);
  fpchl_f.ctr(fpchl_m,2,1);
  hstif_f.ctr(hstif_m,2,1);
  
#define PASS_ADDR(alias,name) (alias) = (name).storage_begin()

  // PAS_ADDR(cartd,dshapex_f);
  cartd = dshapex_f.storage_begin();
  dvolu = wpgdet_v.storage_begin();
  ehist = ehist_f.storage_begin();
  emass = emass_m.storage_begin();
  epmtx = epmtx_m.storage_begin();
  gpcod = gpcod_f.storage_begin();
  props = props_m.storage_begin();
  rmat1 = rmat1_m.storage_begin();
  shape = shape_f.storage_begin();

  PASS_ADDR(stran,stran_f);
  PASS_ADDR(strsg,strsg_f);
  PASS_ADDR(stra0,stra0_f);
  PASS_ADDR(strs0,strs0_f);
  PASS_ADDR(bmatx,bmatx_f);
  PASS_ADDR(desig,desig_m);
  PASS_ADDR(dmatx,dmatx_f);
  PASS_ADDR(dstra,dstra_m);
  PASS_ADDR(presg,presg_m);
  PASS_ADDR(sgtot,sgtot_m);
  PASS_ADDR(sigma,sigma_m);
  PASS_ADDR(tstra,tstra_m);
  PASS_ADDR(tenod,tenod_m);
  PASS_ADDR(dteno,dteno_m);
  PASS_ADDR(pwoel,pwoel_m);
  PASS_ADDR(tgapl,tgapl_m);
  PASS_ADDR(preal,preal_m);
  PASS_ADDR(tenoi,tenoi_m);
  PASS_ADDR(xjacm,xjacm_f);
  PASS_ADDR(bsbar,bsbar_f);
  PASS_ADDR(xjacn,xjacn_f);
  PASS_ADDR(fpchl,fpchl_f);
  PASS_ADDR(hstif,hstif_f);

  thick = &thickness;

  // Connectivities
  lnods = &ICONE(elem,0);

  // ------- ALREADY INVERTED -----------
#define FORTRAN_EXPORT(name) ((FastMat2 &) (name)).storage_begin()

  elcod = FORTRAN_EXPORT(xloc);
  eldis = FORTRAN_EXPORT(state_new);
  bmsig = FORTRAN_EXPORT(res);
  estif = FORTRAN_EXPORT(mat);

  frin30_(cartd,dvolu,ehist,elcod,eldis,emass,epmtx,
	  gpcod,lnods,props,rmat1,shape,stran,strsg,
	  stra0,strs0,thick,
	  bmsig,bmatx,desig,dmatx,dstra,presg,sgtot,
	  sigma,tstra,xjacm,tenod,dteno,
	  pwoel,preal,tgapl,tenoi,bsbar,xjacn,fpchl);

  stif30_(cartd,dvolu,ehist,elcod,eldis,epmtx,gpcod,
	 props,shape,strsg,estif,hstif,emass,
	 bmatx,dmatx,sigma,xjacm,tenod,bsbar,dteno,
	 xjacn);

  // Change sign
  res.scale(-1.);

}

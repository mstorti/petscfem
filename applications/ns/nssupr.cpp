//__INSERT_LICENSE__
/* $Id: nssupr.cpp,v 1.6 2002/09/22 13:18:21 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

#ifdef ROSI_COUPLING_MODULE
#warning Compiling with ROSI_COUPLING_MODULE enabled
extern double AXIAL_ACCELERATION, GLOBAL_TIME;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ns_sup_res::init() {
  int ierr=0;
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,gravity,1.);
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,rho,1.);
  //o Dimension of problem
  TGETOPTDEF_ND(thash,int,ndim,0);
  assert(ndim>0);

  p_indx = ndim+1;
}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)

void ns_sup_res::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & w,FastMat2 & jac) {
  double p, eta;

  // We have (for 3D) first node: {u,v,w,p}
  //                  second node: {eta,lambda,*,*}
  assert(nel==2);
  // Discard k and eps eqs.
  w.set(0.).setel(1.,1,ndim+1,1);
  
#ifdef ROSI_COUPLING_MODULE
//    printf("ns_sup_res: taking t=%f, ax=%f\n",
//  	 GLOBAL_TIME,AXIAL_ACCELERATION);
  double total_axial_acc;
  total_axial_acc = gravity + AXIAL_ACCELERATION;
#define gravity total_axial_acc
#endif

  p = U.get(1,p_indx);
  eta = U.get(2,1);
  r.setel(p-rho*gravity*eta,1);
  jac.setel(1.,1,1,p_indx).setel(-rho*gravity,1,2,1);
}


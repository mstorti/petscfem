//__INSERT_LICENSE__
/* $Id: nssupr.cpp,v 1.1 2001/10/06 23:39:59 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include "rockns.h"

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void ns_sup_res::init() {
  int ierr=0;
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,gravity,1.);
  //o Gravity acceleration. 
  TGETOPTDEF_ND(thash,double,rho,1.);
}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)



void ns_sup_res::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & w,FastMat2 & jac) {
  double p, eta;

  // We have (for 3D) first node: {u,v,w,p}
  //                  second node: {eta,lambda,*,*}
  assert(nel==2);
  // Discard k and eps eqs.
  w.set(0.).setel(1.,1,4,1);
  
  p = U.get(1,4);
  eta = U.get(2,1);
  r.setel(p-eta,1);
  jac.setel(1.,1,4).setel(-1.,2,1);
}


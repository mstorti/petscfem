//__INSERT_LICENSE__
/* $Id: nssupr.cpp,v 1.10 2002/10/13 13:59:46 mstorti Exp $ */

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
double AXIAL_ACCELERATION=DBL_MAX, GLOBAL_TIME;
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

#ifdef ROSI_COUPLING_MODULE
  //o Is being PETSc-FEM called from ROSI?
  TGETOPTDEF_ND(thash,int,called_from_rosi,0);
#endif
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
  double total_axial_acc = gravity;
  // Verify that AXIAL_ACCELERATION has been set (probably in `rosi_hook')
  if (called_from_rosi) {
    assert(AXIAL_ACCELERATION!=DBL_MAX);
    total_axial_acc += AXIAL_ACCELERATION;
    // We could have problems with the free surface if the axal acceleration is too low
    if (total_axial_acc < 0.3 * gravity) total_axial_acc = 0.3 * gravity;
  }
#define gravity total_axial_acc
#endif

  p = U.get(1,p_indx);
  eta = U.get(2,1);
  r.setel(p-rho*gravity*eta,1);
  jac.setel(1.,1,1,p_indx).setel(-rho*gravity,1,2,1);
}


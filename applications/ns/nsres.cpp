//__INSERT_LICENSE__
/* $Id: nsres.cpp,v 1.1 2001/10/16 18:34:21 mstorti Exp $ */

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include <applications/ns/nsi_tet.h>
#include <applications/ns/nssup.h>

extern TextHashTable *GLOBAL_OPTIONS;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void linear_restriction::init() {
  if (!coef) {
    //o The file containing coefficients
    TGETOPTDEF_S(thash,string,coef_file,"coef_file.dat");
    

}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)

void linear_restriction::res(int k,FastMat2 & U,FastMat2 & r,
		       FastMat2 & w,FastMat2 & jac) {
  double p, eta;

  // We have (for 3D) first node: {u,v,w,p}
  //                  second node: {eta,lambda,*,*}
  assert(nel==2);
  // Discard k and eps eqs.
  w.set(0.).setel(1.,1,ndim+1,1);
  
  p = U.get(1,p_indx);
  eta = U.get(2,1);
  r.setel(p-eta,1);
  jac.setel(1.,1,1,p_indx).setel(-1.,1,2,1);
}


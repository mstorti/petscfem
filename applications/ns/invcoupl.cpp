//__INSERT_LICENSE__
//$Id: invcoupl.cpp,v 1.3 2003/02/24 23:19:31 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "./nsi_tet.h"
#include "./adaptor.h"
#include "./invcoupl.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void inviscid_coupling::elemset_init() {
  // NS incompressible
  assert(ndof==ndim+1);
  int ierr;
  //o Dynamic viscosity
  TGETOPTDEF_ND(thash,double,viscosity,1.);

//    tmp1.resize(1,ndim);
//    tmp2.resize(2,nel,ndim);
//    tmp3.resize(3,ndim,nel,ndim);
//    tmp4.resize(4,nel,ndim,nel,ndim);
//    gsopg.resize(2,ndim,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void inviscid_coupling::elemset_end() {
//    tmp1.clear();
//    tmp2.clear();
//    gsopg.clear();
//    tmp3.clear();
//    tmp4.clear();
  tmp.clear();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void inviscid_coupling
::pg_connector(const FastMat2 &xpg,
	       const FastMat2 &state_old_pg,
	       const FastMat2 &grad_state_old_pg,
	       const FastMat2 &state_new_pg,
	       const FastMat2 &grad_state_new_pg,
	       FastMat2 &res_pg,FastMat2 &mat_pg) {

#define tmp1 (tmp(1))
#define tmp2 (tmp(2))
#define tmp3 (tmp(3))
#define tmp4 (tmp(4))
#define gsnpg (tmp(5))
  double fac=-1.;
  // Define residual
  gsnpg.set(grad_state_new_pg);
  gsnpg.is(2,1,ndim);
  // tmp1.prod(normal(),gsnpg,-1,-1,1); // this is wrong!!
  tmp1.prod(gsnpg,normal(),1,-1,-1); // tmp1 is ndim
  tmp2.prod(shape(),tmp1,1,2).scale(-2.0*viscosity*fac);	// tmp2 is nel x ndim
  gsnpg.rs();
  res_pg.set(0.).is(2,1,ndim).set(tmp2).rs();
  // res_pg.print("");

  // Define Jacobian
  tmp3.prod(dshapex(),normal(),1,2,3).scale(2.0*viscosity*fac); // is ndim*nel*ndim
  tmp4.prod(shape(),tmp3,1,2,3,4); // is nel*ndim*nel*ndim
  mat_pg.set(0.).is(2,1,ndim).is(4,1,ndim).set(tmp4).rs();
  // mat_pg.print("");
}

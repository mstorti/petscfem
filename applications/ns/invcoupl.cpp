//__INSERT_LICENSE__
//$Id: invcoupl.cpp,v 1.4 2003/02/27 03:32:41 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/fm2temp.h>

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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void inviscid_coupling::elemset_end() { tmp.clear(); }

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void inviscid_coupling
::pg_connector(const FastMat2 &xpg,
	       const FastMat2 &state_old_pg,
	       const FastMat2 &grad_state_old_pg,
	       const FastMat2 &state_new_pg,
	       const FastMat2 &grad_state_new_pg,
	       FastMat2 &res_pg,FastMat2 &mat_pg) {

#define gsnpg (tmp(0))
  double fac=-1.;
  // Define residual
  gsnpg.set(grad_state_new_pg);
  gsnpg.is(2,1,ndim);
  tmp(1).prod(gsnpg,normal(),1,-1,-1); // tmp1 is ndim
  tmp(2).prod(shape(),tmp(1),1,2).scale(-2.0*viscosity*fac);	// tmp(2) is nel x ndim
  gsnpg.rs();
  res_pg.set(0.).is(2,1,ndim).set(tmp(2)).rs();

  // Define Jacobian
  tmp(3).prod(dshapex(),normal(),1,2,3).scale(2.0*viscosity*fac); // is ndim*nel*ndim
  tmp(4).prod(shape(),tmp(3),1,2,3,4); // is nel*ndim*nel*ndim
  mat_pg.set(0.).is(2,1,ndim).is(4,1,ndim).set(tmp(4)).rs();
}

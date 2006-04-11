//__INSERT_LICENSE__
//$Id: errestim.cpp,v 1.3 2006/04/11 01:31:00 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "errestim.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void error_estimator::
elemset_init() {
  assert(ndof==1);
  int ierr;
  //o Exponent of the norm taken
  TGETOPTDEF_ND(thash,double,norm_expo,2.);
  tmp.resize(1,ndof);
  du.resize(1,ndof);
  G.resize(1,ndof);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void error_estimator::
pg_connector(const FastMat2 &xpg,
	     const FastMat2 &state_old_pg,
	     const FastMat2 &grad_state_old_pg,
	     const FastMat2 &state_new_pg,
	     const FastMat2 &grad_state_new_pg,
	     FastMat2 &res_pg,FastMat2 &mat_pg) {
  assert(ndof==1);
  du.set(state_new_pg).rest(state_old_pg).scale(rec_Dt);
  double g = grad_H.sum_square_all();
  G.setel(g,1);
  tmp.set(G); // G will be the error estimator
  // tmp.print("tmp: ");
  tmp.rest(du);
  res_pg.prod(shape(),tmp,1,2);
  mat_pg.set(0.0);
#if 0
  xpg.print("xpg:");
  H.print("H: ");
  grad_H.print("grad_H: ");
#endif
  
}

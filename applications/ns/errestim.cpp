//__INSERT_LICENSE__
//$Id: errestim.cpp,v 1.1 2006/04/10 22:15:10 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "qharm.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void error_estimator::elemset_init() {
  assert(ndof==1);
  int ierr;
  //o Exponent of the norm taken
  TGETOPTDEF_ND(thash,double,norm_expo,1.);
  //o Specific heat
  TGETOPTDEF_ND(thash,double,rho_Cp,1.);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharm::pg_connector(const FastMat2 &xpg,
			 const FastMat2 &state_old_pg,
			 const FastMat2 &grad_state_old_pg,
			 const FastMat2 &state_new_pg,
			 const FastMat2 &grad_state_new_pg,
			 FastMat2 &res_pg,FastMat2 &mat_pg) {
  assert(ndof==1);
  res_pg
    .prod(dshapex(),grad_state_new_pg,-1,1,-1,2)
    .scale(-conductivity);
  mat_pg
    .ir(2,1)			// Set `dof' indices
    .ir(4,1)
    .prod(dshapex(),dshapex(),-1,1,-1,2) // add prod.
    .scale(conductivity)	// scale by conductivity
    .rs();
  
}

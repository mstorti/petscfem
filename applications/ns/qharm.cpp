//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.4 2002/12/16 04:11:35 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "qharm.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharm::elemset_init() {
  assert(ndof==1);
  int ierr;
  //o Thermal conductivity
  TGETOPTDEF_ND(thash,double,conductivity,1.);
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
  res_pg
    .prod(dshapex,grad_state_new_pg,-1,1,-1,2)
    .scale(-conductivity);
  mat_pg
    .ir(2,1)			// Set `dof' indices
    .ir(4,1)
    .prod(dshapex,dshapex,-1,1,-1,2) // add prod.
    .scale(conductivity)	// scale by conductivity
    .rs();
  
}

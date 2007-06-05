//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "poiboltz.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void poisson_boltzmann::elemset_init() {
  assert(ndof==1);
  int ierr;
  //o 
  TGETOPTDEF_ND(thash,double,ninf,1.0e-6);
  //o 
  TGETOPTDEF_ND(thash,int, z,-1);
  //o 
  TGETOPTDEF_ND(thash,double,eps,8.85e-12);
  //o 
  TGETOPTDEF_ND(thash,double,eps0,80.0);
  //o 
  TGETOPTDEF_ND(thash,double,F,96485.3415);
  //o 
  TGETOPTDEF_ND(thash,double,Tabs,300);
  //o 
  TGETOPTDEF_ND(thash,double,R,8.314472);

 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void poisson_boltzmann::pg_connector(const FastMat2 &xpg,
				     const FastMat2 &state_old_pg,
				     const FastMat2 &grad_state_old_pg,
				     const FastMat2 &state_new_pg,
				     const FastMat2 &grad_state_new_pg,
				     FastMat2 &res_pg,FastMat2 &mat_pg) {
  double psi      = state_new_pg.get(1);
  double sinh_psi = 2*ninf*z*F/(eps*eps0)*sinh(psi*z*F/R/Tabs);
  double cosh_psi = 2*ninf*z*z*F*F/(eps*eps0*R*Tabs)*cosh(psi*z*F/R/Tabs);

  res_pg
    .prod(dshapex(),grad_state_new_pg,-1,1,-1,2)
    .ir(2,1)
    .axpy(shape(), sinh_psi)
    .rs()
    .scale(-1);

  tmp(1)
    .prod(dshapex(),dshapex(),-1,1,-1,2);
    
  tmp(2)
    .prod(shape(),shape(),1,2)
    .scale(cosh_psi);

  mat_pg
    .ir(2,1)
    .ir(4,1)
    .set(tmp(1))
    .add(tmp(2))
    .rs();

  
}

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
  TGETOPTDEF(thash,double,const_1,1.);
  //o 
  TGETOPTDEF(thash,int,   const_2,1);
  //o 
  TGETOPTDEF(thash,double,const_3,1.);
  //o 
  TGETOPTDEF(thash,double,const_4,1.);
  //o 
  TGETOPTDEF(thash,double,const_5,1.);
  //o 
  TGETOPTDEF(thash,double,const_6,1.);

  ninf = const_1;
  z    = const_2;
  eps  = const_3;
  eps0 = const_4; 
  kb   = const_5;
  Tabs = const_6;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void poisson_boltzmann::pg_connector(const FastMat2 &xpg,
				     const FastMat2 &state_old_pg,
				     const FastMat2 &grad_state_old_pg,
				     const FastMat2 &state_new_pg,
				     const FastMat2 &grad_state_new_pg,
				     FastMat2 &res_pg,FastMat2 &mat_pg) {

  double psi      = (double)state_new_pg;
  double sinh_psi = sinh(psi);
  double cosh_psi = cosh(psi);

  res_pg
    .prod(dshapex(),grad_state_new_pg,-1,1,-1,2)
    .scale(-1)
    .axpy(shape(), cosh_psi);

  tmp(1)
    .prod(dshapex(),dshapex(),-1,1,-1,2)
    .scale(1);

  tmp(2)
    .prod(shape(),shape(),1,2)
    .scale(cosh_psi);

  mat_pg
    .ir(2,1)
    .ir(4,1)
    .add(tmp(1)).rest(tmp(2))
    .rs();

  
}

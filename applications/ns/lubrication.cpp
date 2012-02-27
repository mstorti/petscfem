//__INSERT_LICENSE__
//$Id: lubrication.cpp,v 1.12 2007/01/30 19:03:44 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "lubrication.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lubrication::elemset_init() {
  int ierr;
  cond.resize(2,ndof,ndof);
  //o Thermal conductivity
  // read_cond_matrix(thash,"conductivity",ndof,cond);
  cond.eye(1.0);

  C.resize(2,ndof,ndof);
  //o Reaction jacobian matrix
  // read_cond_matrix(thash,"C",ndof,C);
  C.set(0.0);

  Cp.resize(2,ndof,ndof);
  //o Reaction jacobian matrix
  // read_cond_matrix(thash,"Cp",ndof,Cp);
  Cp.set(0.0);

  //o Steady flag
  // TGETOPTDEF(thash,int,steady,0);
  int steady=1;

  //o Time step
  TGETOPTDEF_ND(thash,double,Dt,NAN);
  PETSCFEM_ASSERT0(steady || !isnan(Dt),"Dt is required if not steady");  
  if (steady && isnan(Dt)) Dt=1.0;
  PETSCFEM_ASSERT0(Dt>=0.0,"Dt is required must be non-negative");  

  rec_Dt = (steady? 0.0 : 1.0/Dt);

  //o _T: double[ndof] _N: state_ref _D: null vector 
  // _DOC: Reference state value. _END
  x_ref.resize(1,ndof).set(0.);
  // ierr = get_double(thash,"state_ref",x_ref.storage_begin(),1,ndof);
  // PETSCFEM_ASSERT0(!ierr,"Erro retriveing option state_ref");  

  G.resize(1,ndof).set(1.0);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lubrication::pg_connector(const FastMat2 &xpg,
			 const FastMat2 &state_old_pg,
			 const FastMat2 &grad_state_old_pg,
			 const FastMat2 &state_new_pg,
			 const FastMat2 &grad_state_new_pg,
			 FastMat2 &res_pg,FastMat2 &mat_pg) {
#define tmp1 tmp(1)
#define tmp2 tmp(2)
  tmp1.prod(dshapex(),grad_state_new_pg,-1,1,-1,2);
  res_pg.prod(tmp1,cond,1,-1,-1,2).scale(-1.);
  tmp(6).set(state_new_pg).minus(x_ref);
  tmp(3).prod(C,tmp(6),1,-1,-1);
  tmp(7).prod(shape(),tmp(3),1,2);
  tmp(8).set(state_new_pg).minus(state_old_pg).scale(rec_Dt);
  tmp(9).prod(Cp,tmp(8),1,-1,-1).minus(G);
  tmp(11).prod(shape(),tmp(9),1,2);
  res_pg.minus(tmp(7)).minus(tmp(11));

  tmp2.prod(dshapex(),dshapex(),-1,1,-1,2);
  mat_pg.prod(tmp2,cond,1,3,2,4);
  tmp(4).prod(shape(),shape(),1,2);
  tmp(5).prod(tmp(4),C,1,3,2,4);
  tmp(10).prod(tmp(4),Cp,1,3,2,4);
  mat_pg.add(tmp(5)).axpy(tmp(10),rec_Dt);
}

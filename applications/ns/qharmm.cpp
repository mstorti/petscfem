//__INSERT_LICENSE__
//$Id: qharmm.cpp,v 1.12 2007/01/30 19:03:44 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "qharmm.h"

void read_cond_matrix(TextHashTable *thash, const char *s,
		      int ndof,FastMat2 &cond) {
  vector<double> v;
  const char *line;
  thash->get_entry(s,line);  
  PETSCFEM_ASSERT(line,"option %s is required",s);  
  read_double_array(v,line);
  if (v.size()==1) {
    cond.eye(v[0]);
  } else if (v.size()==(unsigned int)ndof) {
    cond.set(0.).d(2,1).set(&*v.begin()).rs();
  } else if (v.size()==(unsigned int)(ndof*ndof)) {
    cond.set(&*v.begin());
  } else PETSCFEM_ERROR("Number of elements in conductivity line inappropriate\n"
			"entered %d values\n", v.size());  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharmm::elemset_init() {
  int ierr;
  cond.resize(2,ndof,ndof);
  //o Thermal conductivity
  read_cond_matrix(thash,"conductivity",ndof,cond);

  C.resize(2,ndof,ndof);
  //o Reaction jacobian matrix
  read_cond_matrix(thash,"C",ndof,C);

  Cp.resize(2,ndof,ndof);
  //o Reaction jacobian matrix
  read_cond_matrix(thash,"Cp",ndof,Cp);

  //o Steady flag
  TGETOPTDEF(thash,int,steady,0);

  //o Time step
  TGETOPTDEF_ND(thash,double,Dt,NAN);
  PETSCFEM_ASSERT0(steady || !isnan(Dt),"Dt is required if not steady");  
  if (steady && isnan(Dt)) Dt=1.0;
  PETSCFEM_ASSERT0(Dt>=0.0,"Dt is required must be non-negative");  

  rec_Dt = (steady? 0.0 : 1.0/Dt);

  //o _T: double[ndof] _N: state_ref _D: null vector 
  // _DOC: Reference state value. _END
  x_ref.resize(1,ndof).set(0.);
  ierr = get_double(thash,"state_ref",x_ref.storage_begin(),1,ndof);
  PETSCFEM_ASSERT0(!ierr,"Erro retriveing option state_ref");  

  const char *line;
  thash->get_entry("source",line);  
  G.resize(1,ndof).set(0.);
  if (line) {
    vector<double> v;
    read_double_array(v,line);
    if (v.size()==1) {
      G.set(v[0]);
    } else if (v.size()==(unsigned int)ndof) {
      cond.set(&v[0]);
    } else PETSCFEM_ERROR("Number of elements in source line inappropriate\n"
			  "entered %d values\n", v.size()); 
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharmm::pg_connector(const FastMat2 &xpg,
			 const FastMat2 &state_old_pg,
			 const FastMat2 &grad_state_old_pg,
			 const FastMat2 &state_new_pg,
			 const FastMat2 &grad_state_new_pg,
			 FastMat2 &res_pg,FastMat2 &mat_pg) {
#define tmp1 tmp(1)
#define tmp2 tmp(2)
  tmp1.prod(dshapex(),grad_state_new_pg,-1,1,-1,2);
  res_pg.prod(tmp1,cond,1,-1,-1,2).scale(-1.);
  tmp(6).set(state_new_pg).rest(x_ref);
  tmp(3).prod(C,tmp(6),1,-1,-1);
  tmp(7).prod(shape(),tmp(3),1,2);
  tmp(8).set(state_new_pg).rest(state_old_pg).scale(rec_Dt);
  tmp(9).prod(Cp,tmp(8),1,-1,-1).rest(G);
  tmp(11).prod(shape(),tmp(9),1,2);
  res_pg.rest(tmp(7)).rest(tmp(11));

  tmp2.prod(dshapex(),dshapex(),-1,1,-1,2);
  mat_pg.prod(tmp2,cond,1,3,2,4);
  tmp(4).prod(shape(),shape(),1,2);
  tmp(5).prod(tmp(4),C,1,3,2,4);
  tmp(10).prod(tmp(4),Cp,1,3,2,4);
  mat_pg.add(tmp(5)).axpy(tmp(10),rec_Dt);
}

//__INSERT_LICENSE__
//$Id: qharmm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

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
  assert(line);
  read_double_array(v,line);
  if (v.size()==1) {
    cond.eye(v[0]);
  } else if (v.size()==ndof) {
    cond.set(0.).d(2,1).set(v.begin()).rs();
  } else if (v.size()==ndof*ndof) {
    cond.set(v.begin());
  } else PETSCFEM_ERROR("Number of elements in conductivity line inappropriate\n"
			"entered %d values\n", v.size());  
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharmm::elemset_init() {
  int ierr;
  cond.resize(2,ndof,ndof);
  // rho_Cp_m.resize(2,ndof,ndof);
  //o Thermal conductivity
  read_cond_matrix(thash,"conductivity",ndof,cond);
  // read_cond_matrix(thash,"rho_Cp",ndof,rho_Cp_m);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void qharmm::pg_connector(const FastMat2 &xpg,
			 const FastMat2 &state_old_pg,
			 const FastMat2 &grad_state_old_pg,
			 const FastMat2 &state_new_pg,
			 const FastMat2 &grad_state_new_pg,
			 FastMat2 &res_pg,FastMat2 &mat_pg) {
  tmp1.prod(dshapex(),grad_state_new_pg,-1,1,-1,2);
  res_pg.prod(tmp1,cond,1,-1,-1,2).scale(-1.);
  tmp2.prod(dshapex(),dshapex(),-1,1,-1,2);
  mat_pg.prod(tmp2,cond,1,3,2,4);
}

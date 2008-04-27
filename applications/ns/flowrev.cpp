//__INSERT_LICENSE__
//$Id$
#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/util2.h>
#include <src/fastmat2.h>

#include "./nsi_tet.h"
#include "./genload.h"
#include "./flowrev.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_reversal::
start_chunk() {
  int ierr;
  //o Dimension of the problem
  TGETOPTNDEF(thash,int,ndim,none);

  //o Penalization coefficient
  TGETOPTNDEF(thash,double,penal_coef,-1.0);
  PETSCFEM_ASSERT0(penal_coef>=0.0,
                   "Penalization coef must be non-negative");  

  vel_indx=1;
  
  dofs.clear();
  dofs.push_back(ndim+1);

  coefs.clear();
  coefs.push_back(penal_coef);

  refvals.clear();
  refvals.push_back(0.0);
  nor.resize(1,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void flow_reversal::
q(FastMat2 &u,FastMat2 &flux,FastMat2 &jac) {
  nor.set(normal);
  u.is(1,1,ndim);
  tmp.prod(normal,u,-1,-1);
  u.rs();
  double unnor = tmp;
  flux.set(0.0);
  jac.set(0.0);
  double 
    *fluxp = flux.storage_begin(),
    *up = u.storage_begin(),
    *jacp = jac.storage_begin();
  if (unnor<0.0) {
    for (unsigned int j=0; j<dofs.size(); j++) {
      int dof = dofs[j]-1;
      fluxp[dof] = coefs[j]*(up[dof]-refvals[j]);
      jacp[ndof*dof+dof] = coefs[j];
    }
  }
}

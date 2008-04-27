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
  TGETOPTDEF_ND(thash,int,ndim,0);
  PETSCFEM_ASSERT0(ndim>0,"ndim must be positive");  

  //o Dimension of the problem
  TGETOPTDEF_ND(thash,int,vel_indx,1);
  PETSCFEM_ASSERT(vel_indx>=1 && vel_indx+ndim<=ndof,
                  "vel_indx %d out of range, ndof %d",
                  vel_indx,ndof);  
  
  dofs.clear();
  const char *line;
  unsigned int ndofs;
  thash->get_entry("dofs",line);
  if(line) {
    dofs.clear();
    read_int_array(dofs,line);
    ndofs = dofs.size();

    thash->get_entry("coefs",line);
    PETSCFEM_ASSERT0(line,"`coefs' line must be present if `dofs' is");  
    coefs.clear();
    read_double_array(coefs,line);
    PETSCFEM_ASSERT(coefs.size()==ndofs,
                    "`coefs' must have the same length as `dofs'\n"
                    "dofs size %d, coefs size %d",
                    dofs.size(),coefs.size());  
    
    thash->get_entry("refvals",line);
    PETSCFEM_ASSERT0(line,"`refvals' line must be present if `dofs' is");  
    refvals.clear();
    read_double_array(refvals,line);
    PETSCFEM_ASSERT(refvals.size()==ndofs,
                    "`refvals' must have the same length as `dofs'\n"
                    "dofs size %d, refvals size %d",
                    dofs.size(),refvals.size());  
  } else {
    //o Penalization coefficient
    TGETOPTNDEF(thash,double,penal_coef,-1.0);
    PETSCFEM_ASSERT0(penal_coef>=0.0,
                   "Penalization coef must be non-negative");  

    dofs.push_back(ndim+2);

    coefs.clear();
    coefs.push_back(penal_coef);

    refvals.clear();
    refvals.push_back(0.0);
  }
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
      fluxp[dof] = coefs[j]*(refvals[j]-up[dof]);
      jacp[ndof*dof+dof] = coefs[j];
    }
  }
}

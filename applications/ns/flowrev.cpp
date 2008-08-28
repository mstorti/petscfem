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

  //o Field index where the components of the
  //  velocity index start (1-base).
  TGETOPTDEF_ND(thash,int,vel_indx,1);
  PETSCFEM_ASSERT(vel_indx>=1 && vel_indx+ndim<=ndof,
                  "vel_indx %d out of range, ndof %d",
                  vel_indx,ndof);  
  
  //o If this is set, then
  //  the code assumes that a thermal NS problem is being run
  //  i.e. #vel_index=1# and the penalization term is added only
  //  on the temperature field (i.e. #dofs=[ndim+2]#).
  TGETOPTDEF(thash,int,thermal_convection,0);
  
  dofs.clear();
  const char *line;
  unsigned int ndofs;
  //o _T: vector<int>  _N: dofs _D: empty vector
  // _DOC: Field indexes (base-1) for which the flow
  // reversal will be applied _END
  thash->get_entry("dofs",line);
  dofs.clear();
  read_int_array(dofs,line);
  ndofs = dofs.size();

  //o _T: vector<double> _N: coefs _D: empty vector
  // _DOC: Penalization coefficients (artificial film coefficients).
  // The length of #coefs# and #dofs# must be the same. _END
  thash->get_entry("coefs",line);
  PETSCFEM_ASSERT0(line,"`coefs' line must be present if `dofs' is");  
  coefs.clear();
  read_double_array(coefs,line);
  PETSCFEM_ASSERT(coefs.size()==ndofs,
                  "`coefs' must have the same length as `dofs'\n"
                  "dofs size %d, coefs size %d",
                  dofs.size(),coefs.size());  
    
  //o _T: vector<double> _N: refvals _D: empty vector
  // _DOC: The reference values for the unknown.  
  // The length of #coefs# and #dofs# must be the same. _END
  thash->get_entry("refvals",line);
  PETSCFEM_ASSERT0(line,"`refvals' line must be present if `dofs' is");  
  refvals.clear();
  read_double_array(refvals,line);
  if (refvals.size()==0) 
    refvals.resize(ndofs,0);
  else {
    PETSCFEM_ASSERT(refvals.size()==ndofs,
                    "`refvals' must have the same length as `dofs'\n"
                    "dofs size %d, refvals size %d",
                    dofs.size(),refvals.size());  
  }

  if(thermal_convection) {
    PETSCFEM_ASSERT0(dofs.size()==0,
                     "`dofs' option must NOT be set it "
                     "thermal_convection is set.");  
    dofs.push_back(ndim+2);
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

//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "poisson.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void poisson::elemset_init() {
  assert(ndof==1);
  int ierr;
  //o 
  TGETOPTDEF_ND(thash,double,scale_factor,1.0);

  TGETOPTDEF_ND(thash,double,eps,8.85e-12);
  //o 
  TGETOPTDEF_ND(thash,double,eps0,80.0);
  //o 
  TGETOPTDEF_ND(thash,double,F,96485.3415);
  //o 
  tmp.resize(2,nel,ndof);

  int ndim = nodedata->ndim;

  TGETOPTDEF_S_ND(thash,string,concentration_name,concentration);

  bool found; int ncols;

  found = nodedata->get_field(concentration_name,&ncols,&concptr);
  PETSCFEM_ASSERT0(found,"Error");  
  PETSCFEM_ASSERT0(ncols==1,"Error");  
  concol.resize(1,nel);
 

}
void poisson::elem_init() 
{
  int ndim = nodedata->ndim;
  int kel = this->elem;
  int nel = this->nel;

  for (int i=0; i<nel; i++) {
    int n = icone[kel*nel+i] - 1;
     
    double* ci = &concptr[n*1];
    
  
    concol.ir(1,i+1).set(ci).rs();

  }

  
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void poisson::pg_connector(const FastMat2 &xpg,
				     const FastMat2 &state_old_pg,
				     const FastMat2 &grad_state_old_pg,
				     const FastMat2 &state_new_pg,
				     const FastMat2 &grad_state_new_pg,
				     FastMat2 &res_pg,FastMat2 &mat_pg) {

  double concpg = conc.prod(shape(),concol,-1,-1);
    
  tmp.set(concpg).scale(F);
  //double max= tmp.max(1,2);
  //double min= tmp.min(1,2);
  //printf("%f\n ",max);
  //printf("%f\n ",min);
  res_pg
    .prod(dshapex(),grad_state_new_pg,-1,1,-1,2).scale(eps*eps0)
    .add(tmp)
    .scale(-1);


 if (scale_factor != 1.0) {
    if (EVAL_RES) res_pg.scale(scale_factor);
   }
}
int poisson::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			      Dofmap *dofmap,const char *jobinfo,int myrank,
			      int el_start,int el_last,int iter_mode,
			      const TimeData *time_) 
{
  this->nodedata = nodedata;
  return adaptor_pg::assemble(arg_data_v,nodedata,
			      dofmap,jobinfo,myrank,
			      el_start,el_last,iter_mode,
			      time_);
}

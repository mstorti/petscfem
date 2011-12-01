//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "pot_grad.h"
////////---------------------------------------------------------


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void pot_grad::elemset_init() {
  
  int ierr;

  
  int ndim = nodedata->ndim;

  TGETOPTDEF_S_ND(thash,string,potname,potential);

  bool found; int ncols;

  
  found = nodedata->get_field(potname,&ncols,&potptr);
  PETSCFEM_ASSERT0(found,"Error");  
  PETSCFEM_ASSERT0(ncols==1,"Error");  
  potcol.resize(1,nel);
  
}

void pot_grad::elem_init() 
{
  int ndim = nodedata->ndim;
  int kel = this->elem;
  int nel = this->nel;

  for (int i=0; i<nel; i++) {
    int n = icone[kel*nel+i] - 1;
    
    double* p = &potptr[n*1];
    
 
    potcol.ir(1,i+1).set(p).rs();
  }
  
  
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void pot_grad::pg_connector(const FastMat2 &xpg,
				   const FastMat2 &state_old_pg,
				   const FastMat2 &grad_state_old_pg,
				   const FastMat2 &state_new_pg,
				   const FastMat2 &grad_state_new_pg,
				   FastMat2 &res_pg,FastMat2 &mat_pg) {

  int itmp = 1;
#define FM2TMP this->tmp(itmp++)
  
  const FastMat2& ustar = state_new_pg;

  celec.prod(dshapex(),potcol,1,-1,-1).scale(-1);
  
  celec.add(ustar).scale(-1);

 
  if (EVAL_RES) {

   res_pg.prod(shape(), celec, 1, 2);
   
  }

  if (EVAL_MAT) {


    
    FastMat2& massM =FM2TMP;
    massM.prod(shape(),shape(),1,2);
    for (int i=0; i<ndof; i++) {
      mat_pg.ir(2,i+1).ir(4,i+1)
	.set(massM)
	.rs();  
    }
    
  }

}

int pot_grad::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
			      Dofmap *dofmap,const char *jobinfo,int myrank,
			      int el_start,int el_last,int iter_mode,
		       const TimeData *time_) 
		       
{
  this->nodedata = nodedata;
  return adaptor_pg::assemble(arg_data_v,nodedata,
			      dofmap,jobinfo,myrank,
			      el_start,el_last,iter_mode,time_);
}

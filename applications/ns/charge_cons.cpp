//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "charge_cons.h"
////////---------------------------------------------------------
/// implements: div(sigma*grad(fi)+grad(beta)+vel*gamma)=0
/////-------------------<............>-------------------------////////
void charge_cons::elemset_init() {
  int ierr;
  //o 
  int ndim = nodedata->ndim;

  TGETOPTDEF_S_ND(thash,string,sigma_name,conductivity);
  TGETOPTDEF_S_ND(thash,string,beta_name,diffusive_term);
  TGETOPTDEF_S_ND(thash,string,gamma_name,convective_term);
  TGETOPTDEF_S_ND(thash,string,velname,velocity);


  bool found; int ncols;
  
  found = nodedata->get_field(sigma_name,&ncols,&sigmaptr);
  assert(found); assert(ncols==1);
  sigmacol.resize(1,nel);
  
  found = nodedata->get_field(beta_name,&ncols,&betaptr);
  assert(found); assert(ncols==1);
  betacol.resize(1,nel);

  found = nodedata->get_field(gamma_name,&ncols,&gammaptr);
  assert(found); assert(ncols==1);
  gammacol.resize(1,nel);

  found = nodedata->get_field(velname,&ncols,&velptr);
  assert(found); assert(ncols==ndim);
  velcol.resize(2,nel,ndim);
 

  
}

void charge_cons::elem_init() 
{
  int ndim = nodedata->ndim;
  int kel = this->elem;
  int nel = this->nel;

  for (int i=0; i<nel; i++) {
    int n = icone[kel*nel+i] - 1;
    
    double* v = &velptr[n*ndim];
    double* s = &sigmaptr[n*1];
    double* b = &betaptr[n*1];    
    double* g = &gammaptr[n*1];
    
    velcol.ir(1,i+1).set(v).rs();
    sigmacol.ir(1,i+1).set(s).rs();
    betacol.ir(1,i+1).set(b).rs();
    gammacol.ir(1,i+1).set(g).rs();
  }
  
  
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void charge_cons::pg_connector(const FastMat2 &xpg,
				   const FastMat2 &state_old_pg,
				   const FastMat2 &grad_state_old_pg,
				   const FastMat2 &state_new_pg,
				   const FastMat2 &grad_state_new_pg,
				   FastMat2 &res_pg,FastMat2 &mat_pg) {

  int itmp = 1;
#define FM2TMP this->tmp(itmp++)

  double gammapg = gamma.prod(shape(),gammacol,-1,-1);
  
  vel.prod(shape(),velcol,-1,-1,1).scale(-1.0*gammapg);

  grad_beta.prod(dshapex(),betacol,1,-1,-1);

  double sigmapg = sigma.prod(shape(),sigmacol,-1,-1);

  
  

 
  
 
  //////////////////////
 if (EVAL_RES) {
   
  

   

   FastMat2& tmp_res = FM2TMP;

   tmp_res
     .set(grad_state_new_pg).scale(sigmapg)
     .ir(2,1)
     .add(vel)
     .add(grad_beta)
     .rs();


   res_pg.prod(dshapex(),tmp_res,-1,1,-1,2);

}



}

int charge_cons::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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

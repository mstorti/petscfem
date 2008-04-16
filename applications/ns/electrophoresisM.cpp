//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "electrophoresisM.h"
////////---------------------------------------------------------

void read_cond_tensor3(TextHashTable *thash, const char *s,
		      int ndof,FastMat2 &cond) {
  vector<double> v;
  const char *line;
  thash->get_entry(s,line);  
  assert(line);
  read_double_array(v,line);
  if (v.size()==1) {
    cond.eye(v[0]);
  } else if (v.size()==(unsigned int)ndof) {
    cond.set(0.).d(2,1).set(&*v.begin()).rs();
  } else if (v.size()==(unsigned int)(ndof*ndof*ndof)) {
    cond.set(&*v.begin());
  } else PETSCFEM_ERROR("Number of elements in conductivity line inappropriate\n"
			"entered %d values\n", v.size());  
}




//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresisM::elemset_init() {
  //assert(ndof==1);
  int ierr;
  //o 
  TGETOPTDEF_ND(thash,double,R,8.314);
  
  //o 
  TGETOPTDEF_ND(thash,double,T,300);
  //o 
  TGETOPTDEF_ND(thash,double,F,96485.3415);

  diff.resize(2,ndof,ndof);
  //o Diffusion 
  read_cond_matrix(thash,"diffusivity",ndof,diff);
  
  zeff.resize(2,ndof,ndof);
  //o factor electroforetico 
  read_cond_matrix(thash,"mobility",ndof,zeff);

  K.resize(2,ndof,ndof);
  //o Reaccion acoplada 
  read_cond_matrix(thash,"reactions_coef",ndof,K);
  
 //  K1.resize(2,ndof,ndof);
//   //o reaccion simple 
//   read_cond_matrix(thash,"reactions_coef",ndof,K1);

  K1.resize(3,ndof,ndof,ndof);
  //o reaccion externa 
  read_cond_tensor3(thash,"inter_reactions_coef",ndof,K1);

  int ndim = nodedata->ndim;

  TGETOPTDEF_S_ND(thash,string,velname,velocity);
  TGETOPTDEF_S_ND(thash,string,potname,potential);

  bool found; int ncols;

  found = nodedata->get_field(velname,&ncols,&velptr);
  assert(found); assert(ncols==ndim);
  velcol.resize(2,nel,ndim);
  //vel.resize(1,ndim);


  found = nodedata->get_field(potname,&ncols,&potptr);
  assert(found); assert(ncols==1);
  potcol.resize(1,nel);

}

void electrophoresisM::elem_init() 
{
  int ndim = nodedata->ndim;
  int kel = this->elem;
  int nel = this->nel;

  for (int i=0; i<nel; i++) {
    int n = icone[kel*nel+i] - 1;
    
    double* v = &velptr[n*ndim];
    double* p = &potptr[n*1];
    
    velcol.ir(1,i+1).set(v).rs();
    potcol.ir(1,i+1).set(p).rs();
  }
  
  
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresisM::pg_connector(const FastMat2 &xpg,
				   const FastMat2 &state_old_pg,
				   const FastMat2 &grad_state_old_pg,
				   const FastMat2 &state_new_pg,
				   const FastMat2 &grad_state_new_pg,
				   FastMat2 &res_pg,FastMat2 &mat_pg) {

  int itmp = 1;
#define FM2TMP this->tmp(itmp++)
  
  double Dt = glob_param->Dt;
  double alpha = glob_param->alpha;
  double rec_Dt = glob_param->steady ? 0 : (1./Dt);
  

  vel.prod(shape(),velcol,-1,-1,1);

  celec.prod(dshapex(),potcol,1,-1,-1).scale(F/R/T);

 
  const FastMat2& u0 = state_old_pg;
  const FastMat2& u1 = state_new_pg;
  FastMat2& u_star = FM2TMP;
  u_star.set(u1).scale(alpha).axpy(u0,1-alpha);

  const FastMat2& grad_u0 = grad_state_old_pg;
  const FastMat2& grad_u1 = grad_state_new_pg;
  FastMat2& grad_u_star = FM2TMP;
  grad_u_star.set(grad_u1).scale(alpha).axpy(grad_u0,1-alpha);
 
  
  // temporal derivative
  
  FastMat2& du = FM2TMP;
  FastMat2& Du = FM2TMP;
  du.set(u_star).rest(u0).scale(rec_Dt/alpha);
  Du.prod(shape(),du,1,2); 
  
  // convective term
  FastMat2& Cu = FM2TMP;
  FastMat2& tmpCu = FM2TMP;
  tmpCu.prod(vel,grad_u_star,-1,-1,1);
  Cu.prod(shape(),tmpCu,1,2);
  
  //migration term
  FastMat2& Mu = FM2TMP;
  FastMat2& tmpMu0 = FM2TMP;
  FastMat2& tmpMu1 = FM2TMP;
  FastMat2& tmpMu2 = FM2TMP;
  tmpMu0.prod(celec,grad_u_star,-1,-1,1).scale(-1);
  tmpMu1.prod(shape(),tmpMu0,1,2);
  tmpMu2.prod(tmpMu1,zeff,1,-1,-1,2);
  Mu.prod(tmpMu2,diff,1,-1,-1,2);

  
// diffusive term
  
  FastMat2& Lu = FM2TMP;
  FastMat2& tmplu = FM2TMP;
  tmplu.prod(dshapex(),grad_u_star,-1,1,-1,2);
  Lu.prod(tmplu,diff,1,-1,-1,2);

  
 
  
  // single reactive term
  
  FastMat2& Ru = FM2TMP;
  FastMat2& tmpru = FM2TMP;
  tmpru.prod(u_star,K,-1,1,-1);
  Ru.prod(shape(),tmpru,1,2);

  
  
  // conjugate reactive term
  
  FastMat2& Ruc = FM2TMP;
  FastMat2& tmpruc = FM2TMP;
  FastMat2& tmpruc1 = FM2TMP;
  tmpruc.prod(u_star,u_star,1,2);
  tmpruc1.prod(K1,tmpruc,1,-1,-2,-1,-2);
  Ruc.prod(shape(),tmpruc1,1,2);

  // external reactive term 

  


  //////////////////////
  res_pg
    .add(Du)
    .add(Cu)
    .add(Lu)
    .add(Mu)
    .add(Ru)
    .add(Ruc)
    .scale(-1);

  //PetscPrintf(PETSCFEM_COMM_WORLD,"hasta aca llegamos ");

}

int electrophoresisM::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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

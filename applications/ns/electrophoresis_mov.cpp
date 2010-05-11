//__INSERT_LICENSE__
//$Id: qharm.cpp,v 1.5 2003/02/24 00:14:23 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "electrophoresis_mov.h"
////////---------------------------------------------------------
//este elemento permite movilidades distribuidas, sin reaccion
/////////////////////////7




//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresis_mov::elemset_init() {
  //assert(ndof==1);
  int ierr;
  
  TGETOPTDEF_ND(thash,double,supg_fact,1.0);
  TGETOPTDEF_ND(thash,double,m_fact,1.0);
  
 
  diff.resize(2,ndof,ndof);
  read_cond_matrix(thash,"diffusivity",ndof,diff);
  
 
  adveff.resize(2,ndof,ndof);
  read_cond_matrix(thash,"adv_factor",ndof,adveff);

 
  Dm.resize(4,nel,ndof,nel,ndof);
  Cm.resize(4,nel,ndof,nel,ndof);
  Hm.resize(2,ndof,ndof);
  tau_supg.resize(2,ndof,ndof);
  int ndim = nodedata->ndim;
  zeff.resize(2,ndof,ndof);
  zeff.set(0.0);

  ndof_aux.resize(1,ndof);

  TGETOPTDEF_S_ND(thash,string,velname,velocity);
  TGETOPTDEF_S_ND(thash,string,potname,potential);
  TGETOPTDEF_S_ND(thash,string,movname,mobility);
  bool found; int ncols;

  found = nodedata->get_field(velname,&ncols,&velptr);
  assert(found); assert(ncols==ndim);
  velcol.resize(2,nel,ndim);
 


  found = nodedata->get_field(potname,&ncols,&potptr);
  assert(found); assert(ncols==1);
  potcol.resize(1,nel);

  found = nodedata->get_field(movname,&ncols,&movptr);
  assert(found); assert(ncols==ndof);
  movcol.resize(2,nel,ndof);
}

void electrophoresis_mov::elem_init() 
{
  int ndim = nodedata->ndim;
  int kel = this->elem;
  int nel = this->nel;

  for (int i=0; i<nel; i++) {
    int n = icone[kel*nel+i] - 1;
    
    double* v = &velptr[n*ndim];
    double* p = &potptr[n*1];
    double* m = &movptr[n*ndof];

    velcol.ir(1,i+1).set(v).rs();
    potcol.ir(1,i+1).set(p).rs();
    movcol.ir(1,i+1).set(m).rs();
  }
  
  
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresis_mov::pg_connector(const FastMat2 &xpg,
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
  zeff_d.prod(shape(),movcol,-1,-1,1);
 
  for (int ii=1; ii<=ndof; ii++) {
    zeff.setel(zeff_d.get(ii),ii,ii);};
  //  zeff.print();
  //  zeff.eye();
  celec.prod(dshapex(),potcol,1,-1,-1);

 
  const FastMat2& u0 = state_old_pg;
  const FastMat2& u1 = state_new_pg;
  FastMat2& u_star = FM2TMP;
  u_star.set(u1).scale(alpha).axpy(u0,1-alpha);

  const FastMat2& grad_u0 = grad_state_old_pg;
  const FastMat2& grad_u1 = grad_state_new_pg;
  FastMat2& grad_u_star = FM2TMP;
  grad_u_star.set(grad_u1).scale(alpha).axpy(grad_u0,1-alpha);


  ////////calculo de tau_supg
  
  FastMat2& vel_mod    = FM2TMP;
  FastMat2& aux_supg1  = FM2TMP;
  FastMat2& aux_supg2  = FM2TMP;
  FastMat2& aux_supg1a = FM2TMP;
  FastMat2& aux_supg2a = FM2TMP;
  FastMat2& aux_supg4  = FM2TMP;
  FastMat2& aux_supg5  = FM2TMP;
  FastMat2& aux_supg6  = FM2TMP;
  FastMat2& tmp9       = FM2TMP;
  FastMat2& c_supg     = FM2TMP;
  FastMat2& k_supg     = FM2TMP;
  FastMat2& aux_supgc1 = FM2TMP;
  FastMat2& aux_supgk1 = FM2TMP;
  FastMat2& adveff_d   = FM2TMP;

  Hm.eye();
  ndof_aux.set(1.0);

  adveff_d.prod(adveff,ndof_aux,1,-1,-1);
  aux_supg1a.prod(zeff_d,celec,1,2).scale(-1.0);
  aux_supg1.prod(zeff,celec,1,3,2).scale(-1.0);
  aux_supg2a.prod(adveff_d,vel,1,2).add(aux_supg1a);//ndof,ndim
  aux_supg2a.prod(ndof_aux,vel,1,2).add(aux_supg1a);//ndof,ndim
  aux_supg2.prod(adveff,vel,1,3,2).add(aux_supg1);//ndof,ndim,ndof
  aux_supg2.prod(Hm,vel,1,3,2).add(aux_supg1);//ndof,ndim,ndof
  double tol=1.0e-18;
  tau_supg.set(0.0);


  aux_supg4.prod(shape(),aux_supg2,1,2,3,4);//nel,ndof,ndim,ndof
  c_supg.prod(aux_supg4,grad_u_star,1,2,-1,-2,-1,-2);//nel,ndof
  aux_supg5.prod(aux_supg2,dshapex(),1,-1,3,-1,2);//ndof,nel,ndof
  aux_supg6.prod(aux_supg2,grad_u_star,1,-1,-2,-1,-2);//ndof
  k_supg.prod(aux_supg5,aux_supg6,-1,1,2,-1);//nel,ndof

  for (int ii=1; ii<=ndof; ii++) {

    vel_mod.norm_p(aux_supg2a.ir(1,ii), 2); aux_supg2a.rs();
    aux_supgc1.norm_p(c_supg.ir(2,ii),  2); c_supg.rs();
    aux_supgk1.norm_p(k_supg.ir(2,ii),  2); k_supg.rs();

    double tau = 0.0;
    double vel_mod_a = vel_mod;
    FastMat2::branch();
    if(vel_mod_a>tol) {
      FastMat2::choose(0);
      double tau1 = 1./SQ((aux_supgc1+tol*tol)/(aux_supgk1+tol));
      double tau2 = 1./SQ(0.5/(rec_Dt));
      double tau3 = 1./SQ((1./tau1)*SQ(vel_mod)/diff.get(ii,ii));
      tau = 1./sqrt(tau1+tau2+tau3);
    }
    FastMat2::leave();

    tau_supg.setel(tau,ii,ii);
    
  };
  tau_supg.scale(supg_fact);


  //////////////////////
 if (EVAL_RES) {
   //temporal derivative
  FastMat2& du = FM2TMP;
  du.set(u_star).rest(u0).scale(rec_Dt/alpha);
  
// convective term

  FastMat2& Cu = FM2TMP;
  FastMat2& tmpCu = FM2TMP;
  FastMat2& tmpCu1 = FM2TMP;
  tmpCu.prod(vel,u_star,1,2).scale(-1);
  tmpCu1.prod(tmpCu,adveff,1,-1,-1,2);

//migration term

  FastMat2& tmpMu0 = FM2TMP;
  FastMat2& tmpMu1 = FM2TMP;
  tmpMu0.prod(celec,u_star,1,2); 
  tmpMu1.prod(tmpMu0,zeff,1,-1,-1,2);


 // diffusive term
  FastMat2& Lu = FM2TMP;
  FastMat2& tmplu = FM2TMP;
  tmplu.prod(diff,grad_u_star,-1,2,1,-1);
  

  // conjugate reactive term
  

  //supg_term
   FastMat2& E_supg = FM2TMP;
   FastMat2& F_supg = FM2TMP;
   FastMat2& G_supg = FM2TMP;
   FastMat2& H_supg = FM2TMP;

   A_supg.prod(aux_supg2,dshapex(),1,-1,3,-1,2);//ndof,nel,ndof
   // B_supg.prod(K1,u_star,1,2,-1,-1);//ndof,ndof
   //H_supg.set(K);
   //C_supg.prod(H_supg,shape(),1,3,2);//ndof,nel,ndof
   D_supg.set(A_supg);//.rest(C_supg);//ndof,nel,ndof CAMBIAR

   F_supg.prod(aux_supg2,grad_u_star,1,-1,-2,-1,-2)//ndof
     //.add(tmpru)// CAMBIAR
     .add(du);
   E_supg.prod(D_supg,F_supg,2,1,-1,-1);//nel,ndof
   G_supg.prod(tau_supg,E_supg,2,-1,1,-1);//nel,ndof
   

  //res
   FastMat2& Ruc = FM2TMP;
  Ruc.prod(shape(),du,1,2);

  tmplu.add(tmpMu1)
        .add(tmpCu1);
  
  Lu.prod(dshapex(),tmplu,-1,1,-1,2);
 

  res_pg
    .add(Lu)
    .add(Ruc)
    .add(G_supg)
    .scale(-1);

}

if (EVAL_MAT) {

  //temporal derivative
  FastMat2& Dm1j = FM2TMP;
  Dm1j.prod(shape(),shape(),1,2).scale(rec_Dt/alpha);
  Dm.set(0.0);
   for (int ii=1; ii<=ndof; ii++) {
  	    Dm.ir(2,ii).ir(4,ii)
  	      .add(Dm1j).rs();}

// convective term
  FastMat2& tmpCmj = FM2TMP;
  FastMat2& tmp1Cmj = FM2TMP;
  tmpCmj.prod(vel,dshapex(),-1,-1,1).scale(-1);
  tmp1Cmj.prod(shape(),tmpCmj,1,2);

  Cm.prod(tmp1Cmj,adveff,1,3,2,4);


//migration term

  FastMat2& Muj = FM2TMP;
  FastMat2& tmpMu0j = FM2TMP;
  FastMat2& tmpMu1j = FM2TMP;

  
  tmpMu0j.prod(celec,dshapex(),-1,-1,1); 
  
  tmpMu1j.prod(shape(),tmpMu0j,1,2);
  Muj.prod(tmpMu1j,zeff,1,3,2,4);

 // diffusive term
  FastMat2& Luj = FM2TMP;
  FastMat2& tmpluj = FM2TMP;
  tmpluj.prod(dshapex(),diff,1,3,2,4);
  //        .add(tmpMu1j);
  Luj.prod(dshapex(),tmpluj,-1,1,-1,2,3,4);


  // single reactive term

  //  FastMat2& Rucj = FM2TMP;
  //  FastMat2& tmpruj = FM2TMP;
  //  tmpruj.prod(shape(),K,2,1,3);
  //  Rucj.prod(shape(),tmpruj,1,2,3,4);

 // conjugate reactive term
 
 // FastMat2& Rmcj = FM2TMP;
 // FastMat2& tmp1rmcj = FM2TMP;
 // tmp1rmcj.set(B_supg).scale(2.0*alpha/rec_Dt);
 // Rmcj.prod(Dm1j,tmp1rmcj,1,3,2,4);


 //  //  supg term
 
  FastMat2& E_supg = FM2TMP;
  FastMat2& F_supg = FM2TMP;
  FastMat2& G_supg = FM2TMP;
  FastMat2& H_supg = FM2TMP;
  Hm.eye().scale(rec_Dt/alpha);//.add(K);//ndof,ndof CAMBIAR
  
  G_supg.prod(Hm,shape(),1,3,2)//ndof,nel,ndof
    .add(A_supg)//nodf,nel,ndof
    
    ;

  E_supg.prod(D_supg,G_supg,2,1,-1,4,3,-1);//nel,ndof,nel,ndof
  F_supg.prod(E_supg,tau_supg,1,-1,3,4,2,-1);
 
  //////////
  mat_pg
    .add(Dm)
    .add(Cm)
    // .add(Rucj)
    .add(Luj)
    .add(Muj)
    //  .add(Rmcj)
    .add(F_supg)
    ;
    
}
  /* res & mat scaling */
  if (m_fact != 1.0) {
    if (EVAL_RES) res_pg.scale(m_fact);
    if (EVAL_MAT) mat_pg.scale(m_fact);
  }
  



}

int electrophoresis_mov::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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


// dimensiones de los vectores y matrices:
// shape():nel
// dshapex():ndim,nel
// celec:ndim
// u_star:ndof
// grad_ustar:ndim,ndof
// K, zeff,diff:ndof,ndof
// K1: ndof,ndof,ndof
// res_pg :nel,ndof
// mat_pg: nel, ndof,nel,ndof

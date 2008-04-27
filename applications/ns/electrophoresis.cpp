//__INSERT_LICENSE__
//$Id$

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "electrophoresis.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresis::elemset_init() {
  //assert(ndof==1);
  int ierr;
  //o 
  TGETOPTDEF_ND(thash,double,N,6.022e23);
  //o 
  TGETOPTDEF_ND(thash,double,el,1.60217646e-19);
  //o 
  TGETOPTDEF_ND(thash,double,mu,1e-3);
  //no son necesariamente iguales para todas las especies
  //o 
  TGETOPTDEF_ND(thash,double,rm,1e-9);
 
  //o 
  TGETOPTDEF_ND(thash,double,diff,1.0);
  //o 
  TGETOPTDEF_ND(thash,double,r,1.0);
  //o 
  TGETOPTDEF_ND(thash,int, z,-1);


  int ndim = nodedata->ndim;

  TGETOPTDEF_S_ND(thash,string,velname,velocity);
  TGETOPTDEF_S_ND(thash,string,potname,potential);

  bool found; int ncols;

  found = nodedata->get_field(velname,&ncols,&velptr);
  assert(found); assert(ncols==ndim);
  velcol.resize(2,nel,ndim);
  vel.resize(1,ndim);


  found = nodedata->get_field(potname,&ncols,&potptr);
  assert(found); assert(ncols==1);
  //potcol.resize(2,nel,1);
  potcol.resize(1,nel);

}

void electrophoresis::elem_init() 
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

  //velcol.print("velcol");
  //potcol.print("potcol");
}


//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void electrophoresis::pg_connector(const FastMat2 &xpg,
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
  celec.prod(dshapex(),potcol,1,-1,-1).scale(-z*el/6.0/3.141592/mu/rm);
  vel.add(celec);
    
  const FastMat2& u0 = state_old_pg;
  const FastMat2& u1 = state_new_pg;
  FastMat2& u_star = FM2TMP;
  u_star.set(u1).scale(alpha).axpy(u0,1-alpha);

  const FastMat2& grad_u0 = grad_state_old_pg;
  const FastMat2& grad_u1 = grad_state_new_pg;
  FastMat2& grad_u_star = FM2TMP;
  grad_u_star.set(grad_u1).scale(alpha).axpy(grad_u0,1-alpha);
  
  // material derivative
  FastMat2& du = FM2TMP;
  FastMat2& dmatu = FM2TMP;
  FastMat2& Du = FM2TMP;
  du.set(u_star).rest(u0);
  dmatu
    .prod(vel,grad_u_star,-1,-1,1)
    .axpy(du, rec_Dt/alpha);
  Du.prod(shape(),dmatu,1,2);

  // diffusive term
  FastMat2& Lu = FM2TMP;
  Lu.prod(dshapex(),grad_u_star,-1,1,-1,2).scale(diff);
 
  // reactive term
  FastMat2& Ru = FM2TMP;
  Ru.prod(shape(),u_star,1,2).scale(r);

  res_pg
    .add(Du)
    .add(Lu)
    .add(Ru)
    .scale(-1);


#if 0
  tmp(1)
    .prod(dshapex(),dshapex(),-1,1,-1,2);
    
  tmp(2)
    .prod(shape(),shape(),1,2)
    .scale(cosh_psi);

  mat_pg
    .ir(2,1)
    .ir(4,1)
    .set(tmp(1))
    .add(tmp(2))
    .rs();
#endif
}

int electrophoresis::assemble(arg_data_list &arg_data_v,Nodedata *nodedata,
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

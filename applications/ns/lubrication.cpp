//__INSERT_LICENSE__
//$Id: lubrication.cpp,v 1.12 2007/01/30 19:03:44 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/util2.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "lubrication.h"

lubrication *lubrication_p=NULL;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lubrication::lub_init() {
  if (!lubrication_p) {
    lubrication_p = this;
    int ierr;
    PETSCFEM_ASSERT0(ndof==1,"lubrication requires...");  
    PETSCFEM_ASSERT0(ndim==2,"lubrication requires...");  

    cond.resize(2,ndof,ndof);
    //o Thermal conductivity
    // read_cond_matrix(thash,"conductivity",ndof,cond);
    cond.eye(1.0);

    C.resize(2,ndof,ndof);
    //o Reaction jacobian matrix
    // read_cond_matrix(thash,"C",ndof,C);
    C.set(0.0);

    Cp.resize(2,ndof,ndof);
    //o Reaction jacobian matrix
    // read_cond_matrix(thash,"Cp",ndof,Cp);
    Cp.set(0.0);

    //o Steady flag
    // TGETOPTDEF(thash,int,steady,0);
    int steady=1;

    //o Time step
    TGETOPTDEF_ND(thash,double,Dt,NAN);
    PETSCFEM_ASSERT0(steady || !isnan(Dt),"Dt is required if not steady");  
    if (steady && isnan(Dt)) Dt=1.0;
    PETSCFEM_ASSERT0(Dt>=0.0,"Dt is required must be non-negative");  

    rec_Dt = (steady? 0.0 : 1.0/Dt);

    //o _T: double[ndof] _N: state_ref _D: null vector 
    // _DOC: Reference state value. _END
    x_ref.resize(1,ndof).set(0.);
    // ierr = get_double(thash,"state_ref",x_ref.storage_begin(),1,ndof);
    // PETSCFEM_ASSERT0(!ierr,"Erro retriveing option state_ref");  

    G.resize(1,ndof).set(1.0);

    TGETOPTDEF_ND(thash,int,use_cavitation_model,0);

#define PFDBREQ(name)                                           \
    TGETOPTDEF_ND(thash,double,name,NAN);                       \
    PETSCFEM_ASSERT0(!isnan(name),#name " is required!!");  
  
    PFDBREQ(rho);
    PFDBREQ(viscosity);
    nu = viscosity/rho;
    PFDBREQ(Omega0);
    PFDBREQ(Omega1);

    PFDBREQ(L);
    PFDBREQ(R);
    PFDBREQ(c);
    double epsil;
    PFDBREQ(epsil);
  
    TGETOPTDEF_ND(thash,double,e0x,0.0);
    TGETOPTDEF_ND(thash,double,e0y,0.0);
    TGETOPTDEF_ND(thash,double,e1x,0.0);
    TGETOPTDEF_ND(thash,double,e1y,0.0);

    TGETOPTDEF_ND(thash,double,e0xdot,0.0);
    TGETOPTDEF_ND(thash,double,e0ydot,0.0);
    TGETOPTDEF_ND(thash,double,e1xdot,0.0);
    TGETOPTDEF_ND(thash,double,e1ydot,0.0);

  } else {
    PETSCFEM_ASSERT0(lubrication_p==this,
                     "Currently there can be only an instance of lubrication.");  
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void lubrication::elemset_init() { 
  lub_init();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void lubrication::pg_connector(const FastMat2 &xpg,
			 const FastMat2 &state_old_pg,
			 const FastMat2 &grad_state_old_pg,
			 const FastMat2 &state_new_pg,
			 const FastMat2 &grad_state_new_pg,
			 FastMat2 &res_pg,FastMat2 &mat_pg) {
  double rhs,
    xx = xpg.get(1),
    phi = xx/R,
    dex = e1x-e0x,
    dey = e1y-e0y,
    h = c-dex*cos(phi)-dey*sin(phi),
    kond = rho*CB(h)/(12.0*viscosity),
    Omega = (Omega1+Omega0)/2.0,
    dexdot = e1xdot-e0xdot,
    deydot = e1ydot-e0ydot;
  rhs = Omega*(-dex*sin(phi)+dey*cos(phi))
    + (dexdot*cos(phi)+deydot*sin(phi));
  rhs *= rho;
  cond.eye(kond);
  G.set(rhs);

#define tmp1 tmp(1)
#define tmp2 tmp(2)
  tmp1.prod(dshapex(),grad_state_new_pg,-1,1,-1,2);
  res_pg.prod(tmp1,cond,1,-1,-1,2).scale(-1.);
  tmp(6).set(state_new_pg).minus(x_ref);
  tmp(3).prod(C,tmp(6),1,-1,-1);
  tmp(7).prod(shape(),tmp(3),1,2);
  tmp(8).set(state_new_pg).minus(state_old_pg).scale(rec_Dt);
  tmp(9).prod(Cp,tmp(8),1,-1,-1).minus(G);
  tmp(11).prod(shape(),tmp(9),1,2);
  res_pg.minus(tmp(7)).minus(tmp(11));

  tmp2.prod(dshapex(),dshapex(),-1,1,-1,2);
  mat_pg.prod(tmp2,cond,1,3,2,4);
  tmp(4).prod(shape(),shape(),1,2);
  tmp(5).prod(tmp(4),C,1,3,2,4);
  tmp(10).prod(tmp(4),Cp,1,3,2,4);
  mat_pg.add(tmp(5)).axpy(tmp(10),rec_Dt);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void lub_force_integrator::init() { 
  PETSCFEM_ASSERT0(lubrication_p,
                   "lubrication elemeset instance should be initialized already");  
  assert(gather_length>=2); 
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void lub_force_integrator::set_pg_values(vector<double> &pg_values,FastMat2 &u,
                   FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
                   double wpgdet,double time) {
  lubrication_p->set_pg_values(pg_values,u,uold,xpg,Jaco,wpgdet,time);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void lubrication::set_pg_values(vector<double> &pg_values,FastMat2 &u,
                                FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
                                double wpgdet,double time) {
  double 
    xx = xpg.get(1),
    P = u.get(1),
    phi = xx/R;
  if (use_cavitation_model && P<0.0) P = 0.0;
  // Just the area so far...
  pg_values[0] = -wpgdet*P*cos(phi);
  pg_values[1] = -wpgdet*P*sin(phi);
}

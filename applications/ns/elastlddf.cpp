//__INSERT_LICENSE__
//$Id: elastld.cpp,v 1.17 2007/01/30 19:03:44 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "nsgath.h"
#include "elastlddf.h"

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define MAXPROPS 100

#define USE_YOUNG_PER_ELEM

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_df::init() {

  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);

  int ierr, iprop=0;
  Young_modulus_indx = iprop; 

#ifdef USE_YOUNG_PER_ELEM
  //o Young modulus
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "Young_modulus",1);
  nprops = iprop;

  //o If Young modulus is entered by element, then this 
  //  factor affects the spatial dependence entered in the 
  //  per-element table. 
  TGETOPTDEF_ND(thash,double,Young_modulus_fac,1.0);
  assert(Young_modulus_fac>0.);
#else
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);
#endif

  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);

  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;

  //o Damping coefficient
  TGETOPTDEF_ND(thash,double,cdamp,0.);

  //o Use new formulation (swap eqs, and rewrite acceleration)
  TGETOPTDEF_ND(thash,int,use_new_form,1);

  G_body.resize(1,ndim).set(0.);
  const char *line;
  vector<double> G_body_v;
  thash->get_entry("G_body",line);
  if(line) {
    read_double_array(G_body_v,line);
    assert(G_body_v.size()==(unsigned int)ndim);
    G_body.set(&G_body_v[0]);
  }

  // Dos opciones para imprimir
  // printf("rec_Dt: %d\n",rec_Dt);
  // SHV(rec_Dt);
  assert(!(rec_Dt>0. && rho==0.));
  assert(ndof==ndim);

  ntens = ndim*(ndim+1)/2;
  nen = nel*ndim;
  
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  grad_u.resize(2,ndim,ndim);
  F.resize(2,ndim,ndim);
  ustar.resize(2,nel,ndim);
  Id.resize(2,ndim,ndim).eye();

  xmh.resize(2,nel,ndim);
  xph.resize(2,nel,ndim);

  // I don't know right if this will work for
  // ndim==3 (may be yes)
  // assert(ndim==2);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_df::
before_chunk(const char *jobinfo) {
  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_mat_res);

  PETSCFEM_ASSERT(comp_prof || comp_mat_res,
                  "Only jobinfo=\"comp_mat_res\" processed. \n"
                  "Received unrecognized jobinfo %s\n",jobinfo);  

  int ierr;
  if (comp_mat_res) {
    state_mh_argh = get_arg_handle("state_mh",
                                   "No handle for `state_mh'\n");
    state_ph_argh = get_arg_handle("state_ph",
                                   "No handle for `state_ph'\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_df
::element_connector(const FastMat2 &xloc,
                    const FastMat2 &state_old,
                    const FastMat2 &state_new,
                    FastMat2 &res,FastMat2 &mat) {
  res.set(0.);
  mat.set(0.);
  get_vals(state_mh_argh,xmh);
  get_vals(state_ph_argh,xph);
  FMSHV(xmh);
  FMSHV(xph);
  PetscFinalize();
  exit(0);
    
#ifdef USE_YOUNG_PER_ELEM
  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));
  double Ex = *(propel.buff()+Young_modulus_indx);
  E = Ex*Young_modulus_fac;
#if 0
  if (rand()%1000==0) {
    xloc.print("xloc:");
    printf("Ex %f, Young_modulus_fac %f\n",Ex,Young_modulus_fac);
  }
#endif
#endif

  // printf("element %d, Young %f\n",elem,E);

  lambda = nu*E/((1+nu)*(1-2*nu));
  mu = E/2/(1+nu);

  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
    dshapex.rs();
    res.rs();

    shape.ir(2,ipg+1);
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco = Jaco.det();
    if (detJaco <= 0.) {
      detj_error(detJaco,elem);
      set_error(1);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);
    iJaco.inv(Jaco);
    dshapex.prod(iJaco,dshapexi,1,-1,-1,2);

    tmp3.set(state_new);
    tmp3.is(2,1,ndim);
    xnew.set(tmp3);
    tmp3.rs().is(2,ndim+1,2*ndim);
    vnew.set(tmp3);
    tmp3.rs();

    tmp3.set(state_old);
    tmp3.is(2,1,ndim);
    xold.set(tmp3);
    tmp3.rs().is(2,ndim+1,2*ndim);
    vold.set(tmp3);
    tmp3.rs();

    ustar.set(xnew).scale(alpha).axpy(xold,1-alpha);
    vstar.set(vnew).scale(alpha).axpy(vold,1-alpha);

    grad_u.prod(ustar,dshapex,-1,1,2,-1);
    F.set(Id).add(grad_u);
    strain.prod(F,F,-1,1,-1,2).rest(Id).scale(0.5);
    tmp5.ctr(strain,-1,-1);
    double trE = tmp5;
    tmp4.set(Id).scale(trE*lambda).axpy(strain,2*mu);
    stress.prod(F,tmp4,1,-1,-1,2);

    // Velocity from states dxdt = (xnew-xold)/dt
    dxdt.set(xnew).rest(xold).scale(rec_Dt);

    // Inertia term
    if (use_new_form) {
      // In this form the residual is [Rmom; Rvel] and
      // acceleration is computed from displacements.  It is
      // better conditioned (I guess), the resulting Jacobian
      // is [I/Dt^2+K,0; -I/Dt,I]
      vnew1.set(dxdt).axpy(vold,-(1.0-alpha)).scale(1.0/alpha);
      a.set(vnew1).rest(vold).scale(rec_Dt)
        .axpy(vstar,cdamp);
    } else {
      // In this form the residual is [Rvel; Rmom]
      // and acceleration is computed from velocities
      // only.  It is bad conditioned, the resulting Jacobian
      // is [I/Dt,-I; I/Dt, K]
      a.set(vnew).rest(vold).scale(rec_Dt)
        .axpy(vstar,cdamp);
    }
    tmp.prod(shape,a,-1,-1,1).rest(G_body);
    tmp2.prod(shape,tmp,1,2);
    res.is(2,ndim+1,2*ndim).axpy(tmp2,-wpgdet*rho);

    // Elastic force residual computation
#if 1
    res_pg.prod(dshapex,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet).rs();
#else
    tmp6.prod(ustar,shape,-1,1,-1);
    res_pg.prod(shape,tmp6,1,2);
    res.axpy(res_pg,1.0).rs();
#endif
    
    // Eqs. for displacements: (xnew-xold)/dt - vstar = 0
    dv.set(dxdt).rest(vstar);
    tmp.prod(shape,dv,-1,-1,1);
    tmp2.prod(shape,tmp,1,2);
    res.is(2,1,ndim).axpy(tmp2,-wpgdet).rs();

  }
  shape.rs();
  res.rs();

  if (use_new_form) {
    // Swap residual components
    // [Rvel; Rmom] -> [Rmom; Rvel]
    res.is(2,1,ndim);
    tmp7.set(res);
    res.rs().is(2,ndim+1,2*ndim);
    tmp8.set(res);
    res.set(tmp7).scale(-1.0);
    res.rs().is(2,1,ndim).set(tmp8).rs();
  }
    
}

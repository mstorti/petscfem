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
#include "elastld.h"

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)
#define MAXPROPS 100

#define USE_YOUNG_PER_ELEM

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity::init() {

  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);

  int ierr;
  //o Use displacement only based formulation. This is normally
  //  set in the struct_main() routine, NOT by the user in the
  //  user data file. 
  TGETOPTDEF_ND(thash,int,use_displacement_formulation,0);

  newmark_gamma = newmark_beta = NAN;
  if (use_displacement_formulation) {
    //o Newmark gamma parameter
    TGETOPTDEF_ND(thash,double,newmark_gamma,0.5);

    //o Newmark beta parameter
    TGETOPTDEF_ND(thash,double,newmark_beta,0.25);
  }

  int iprop=0;
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
  PETSCFEM_ASSERT(!(rec_Dt>0. && rho==0.),
                  "If not steady, then rho must be positive. rec_Dt %f, rho %f",
                  rec_Dt,rho);
  if (!use_displacement_formulation) {
    PETSCFEM_ASSERT(ndof==2*ndim,
                    "If the displacement-velocity formulation is used, then ndof==2*ndim."
                    "Entered ndof %d, ndim %d",ndof,ndim);  
  } else {
    PETSCFEM_ASSERT(ndof==ndim,
                    "If the displacement only formulation is used, then ndof==ndim."
                    "Entered ndof %d, ndim %d",ndof,ndim);  
  }

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

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity::
before_chunk(const char *jobinfo) {
  GET_JOBINFO_FLAG(comp_prof);
  GET_JOBINFO_FLAG(comp_mat_res);

  PETSCFEM_ASSERT(comp_prof || comp_mat_res,
                  "Only jobinfo=\"comp_mat_res\" processed. \n"
                  "Received unrecognized jobinfo %s\n",jobinfo);  

  int ierr;
  if (comp_mat_res && use_displacement_formulation) {
    res_h = get_arg_handle("res","No handle for `res'\n");
    mat_h = get_arg_handle("A","No handle for `A'\n");
    state_mh_argh = get_arg_handle("state_mh",
                                   "No handle for `state_mh'\n");
    state_ph_argh = get_arg_handle("state_ph",
                                   "No handle for `state_ph'\n");
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity
::element_connector(const FastMat2 &xloc,
                    const FastMat2 &state_old,
                    const FastMat2 &state_new,
                    FastMat2 &res,FastMat2 &mat) {

  if (use_displacement_formulation) {
    element_connector_df(xloc,state_old,state_new,res,mat);
    return;
  }

  res.set(0.);
  mat.set(0.);

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
    strain.prod(F,F,-1,1,-1,2).minus(Id).scale(0.5);
    tmp5.ctr(strain,-1,-1);
    double trE = tmp5;
    tmp4.set(Id).scale(trE*lambda).axpy(strain,2*mu);
    stress.prod(F,tmp4,1,-1,-1,2);

    // Velocity from states dxdt = (xnew-xold)/dt
    dxdt.set(xnew).minus(xold).scale(rec_Dt);

    // Inertia term
    if (use_new_form) {
      // In this form the residual is [Rmom; Rvel] and
      // acceleration is computed from displacements.  It is
      // better conditioned (I guess), the resulting Jacobian
      // is [I/Dt^2+K,0; -I/Dt,I]
      vnew1.set(dxdt).axpy(vold,-(1.0-alpha)).scale(1.0/alpha);
      a.set(vnew1).minus(vold).scale(rec_Dt)
        .axpy(vstar,cdamp);
    } else {
      // In this form the residual is [Rvel; Rmom]
      // and acceleration is computed from velocities
      // only.  It is bad conditioned, the resulting Jacobian
      // is [I/Dt,-I; I/Dt, K]
      a.set(vnew).minus(vold).scale(rec_Dt)
        .axpy(vstar,cdamp);
    }
    tmp.prod(shape,a,-1,-1,1).minus(G_body);
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
    dv.set(dxdt).minus(vstar);
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_load
::init() {
  
  int ierr;
  //o Use the displacement only formulation
  TGETOPTDEF(thash,int,use_displacement_formulation,0);

#define MAXPROPS 100
  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);

  force.resize(1,ndim).set(0.);
  const char *line2;
  vector<double> force_v;
  thash->get_entry("Tractions",line2);
  if(line2) {
    read_double_array(force_v,line2);
    assert(force_v.size()==(unsigned int)ndim);
    force.set(&force_v[0]);
    //printf("Readed Tractions Fx=%f Fy=%f Fz=%f  \n",force_v[0],force_v[1],force_v[2]);
  }
  
  int iprop=0;
  pressure_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "pressure",1);
  nprops = iprop;

  //o Poisson ratio
  TGETOPTDEF(thash,double,defo_fac,1.);
  //o Use new formulation (swap eqs, and rewrite acceleration)
  TGETOPTDEF_ND(thash,int,use_new_form,1);
  
  nor.resize(1,ndim);
  Jaco.resize(2,ndimel,ndim);
  tmp.resize(2,nel,ndim);
  tmp1.resize(2,nel,ndim);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_load
::element_connector(const FastMat2 &xloc,
		    const FastMat2 &state_old,
		    const FastMat2 &state_new,
		    FastMat2 &res,FastMat2 &mat) {

  res.set(0.);
  mat.set(0.);
  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));

  double pressure = *(propel.buff()+pressure_indx);
  
  if (use_displacement_formulation) {
    xstar.set(xloc).axpy(state_new,defo_fac);
  } else {
    if (use_new_form) res.is(2,1,ndim);
    else res.is(2,ndim+1,2*ndim);

    xstar.set(xloc);

    state.set(state_old);
    state.is(2,1,ndim);
    xstar.axpy(state,defo_fac*(1-alpha));
    state.rs();

    state.set(state_new);
    state.is(2,1,ndim);
    xstar.axpy(state,defo_fac*alpha);
    state.rs();
  }

  for (int ipg=0; ipg<npg; ipg++) {

    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    shape.ir(2,ipg+1);
    Jaco.prod(dshapexi,xstar,1,-1,-1,2);
    double detJaco = Jaco.detsur(&nor);
    if (detJaco<=0.) {
      detj_error(detJaco,elem);
      set_error(1);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);

    tmp.prod(shape,nor,1,2).scale(-wpg.get(ipg+1)*pressure);
    tmp1.prod(shape,force,1,2).scale(wpgdet);
    res.add(tmp).add(tmp1);
  }
  dshapexi.rs();
  shape.rs();
  res.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity
::element_connector_df(const FastMat2 &xloc,
                       const FastMat2 &state_old,
                       const FastMat2 &state_new,
                       FastMat2 &res,FastMat2 &mat) {
  res.set(0.);
  mat.set(0.);
  xnew.set(state_new);
  xold.set(state_old);
  get_vals(state_mh_argh,xmh);
  get_vals(state_ph_argh,xph);
    
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

  double Dt = 1.0/rec_Dt, Dt2 = Dt*Dt;

  vold.set(xph).minus(xmh).scale(rec_Dt);
  aold.set(xph).axpy(xold,-2.0).add(xmh).scale(4.0*rec_Dt*rec_Dt);
  anew.set(xnew).minus(xold).axpy(vold,-Dt).scale(2.0/Dt2)
    .axpy(aold,-(1.0-2.0*newmark_beta)).scale(1.0/(2.0*newmark_beta));
  vnew.set(vold)
    .axpy(aold,(1.0-newmark_gamma)*Dt)
    .axpy(anew,newmark_gamma*Dt);

#if 0
  if (prtb_index()==0) {
    printf("------------------\nrec_Dt %f\n",rec_Dt);
    FMSHV(xmh);
    FMSHV(xold);
    FMSHV(xph);
    FMSHV(xnew);
    FMSHV(vold);
    FMSHV(aold);
    FMSHV(vnew);
    FMSHV(anew);
    printf("-----------------------------------\n");
  }
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


    grad_u.prod(xnew,dshapex,-1,1,2,-1);
    F.set(Id).add(grad_u);
    strain.prod(F,F,-1,1,-1,2).minus(Id).scale(0.5);
    tmp5.ctr(strain,-1,-1);
    double trE = tmp5;
    tmp4.set(Id).scale(trE*lambda).axpy(strain,2*mu);
    stress.prod(F,tmp4,1,-1,-1,2);

    tmp.prod(shape,anew,-1,-1,1).minus(G_body);
    tmp2.prod(shape,tmp,1,2);
    res.axpy(tmp2,-wpgdet*rho);

    // Elastic force residual computation
    res_pg.prod(dshapex,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet).rs();
    
  }
  shape.rs();
  res.rs();
  // export_vals(res_h,res);
  // export_vals(mat_h,mat);
}

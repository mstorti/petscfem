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

  //o Newmark gamma parameter
  TGETOPTDEF_ND(thash,double,newmark_gamma,0.5);

  //o Newmark beta parameter
  TGETOPTDEF_ND(thash,double,newmark_beta,0.25);

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
  PETSCFEM_ASSERT0(rec_Dt>0.0,"Dt must be positive and finite.");

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
                  "Only jobinfos=\"comp_mat_res\" or \"comp_prof\" processed. \n"
                  "Received unrecognized jobinfo %s\n",jobinfo);  

  int ierr;
  if (comp_mat_res) {
    res_h = get_arg_handle("res","No handle for `res'\n");
    mat_h = get_arg_handle("A","No handle for `A'\n");
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

  vold.set(xph).rest(xmh).scale(rec_Dt);
  aold.set(xph).axpy(xold,-2.0).add(xmh).scale(rec_Dt*rec_Dt);
  anew.set(xnew).rest(xold).axpy(vold,-Dt).scale(2.0/Dt2)
    .axpy(aold,1.0-2.0*newmark_beta).scale(1.0/(2.0*newmark_beta));
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
    strain.prod(F,F,-1,1,-1,2).rest(Id).scale(0.5);
    tmp5.ctr(strain,-1,-1);
    double trE = tmp5;
    tmp4.set(Id).scale(trE*lambda).axpy(strain,2*mu);
    stress.prod(F,tmp4,1,-1,-1,2);

    tmp.prod(shape,anew,-1,-1,1).rest(G_body);
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

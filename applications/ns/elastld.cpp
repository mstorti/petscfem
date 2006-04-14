//__INSERT_LICENSE__
//$Id: elastld.cpp,v 1.12 2006/04/14 13:44:51 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "nsgath.h"
#include "elastld.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity::init() {

  int ierr;
  //o Young modulus
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);

  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);

  lambda = nu*E/((1+nu)*(1-2*nu));
  mu = E/2/(1+nu);

  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;

  //o Damping coefficient
  TGETOPTDEF_ND(thash,double,cdamp,0.);

  G_body.resize(1,ndim).set(0.);
  const char *line;
  vector<double> G_body_v;
  thash->get_entry("G_body",line);
  if(line) {
    read_double_array(G_body_v,line);
    assert(G_body_v.size()==ndim);
    G_body.set(&G_body_v[0]);
  }

  // Dos opciones para imprimir
  // printf("rec_Dt: %d\n",rec_Dt);
  // SHV(rec_Dt);
  assert(!(rec_Dt>0. && rho==0.));
  assert(ndof==2*ndim);

  ntens = ndim*(ndim+1)/2;
  nen = nel*ndim;
  
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  mass_pg.resize(2,nel,nel);
  grad_u.resize(2,ndim,ndim);
  F.resize(2,ndim,ndim);
  ustar.resize(2,nel,ndim);
  Id.resize(2,ndim,ndim).eye();

  // I don't know right if this will work for
  // ndim==3 (may be yes)
  assert(ndim==2);

}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  res.set(0.);
  mat.set(0.);

  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
    dshapex.rs();
    res.rs();

    shape.ir(2,ipg+1);
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco = Jaco.det();
    if (detJaco <= 0.) {
      printf("Jacobian of element %d is negative or null\n"
	     " Jacobian: %f\n",elem,detJaco);
      PetscFinalize();
      exit(0);
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

    // Inertia term
    a.set(vnew).rest(vold).scale(rec_Dt)
      .axpy(vstar,cdamp);
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
    
    mass_pg.prod(shape,shape,1,2).scale(wpgdet*rec_Dt*rho/alpha);
    for (int k=1; k<=ndim; k++) 
      mat.ir(2,ndim+k).ir(4,ndim+k).add(mass_pg);
    mat.rs();

    // Eqs. for displacements: (xnew-xold)/dt - vstar = 0
    dv.set(xnew).rest(xold).scale(rec_Dt).rest(vstar);
    tmp.prod(shape,dv,-1,-1,1);
    tmp2.prod(shape,tmp,1,2);
    res.is(2,1,ndim).axpy(tmp2,-wpgdet);

    mass_pg.prod(shape,shape,1,2).scale(wpgdet);

    for (int k=1; k<=ndim; k++) {
      mat.ir(2,k).ir(4,k).axpy(mass_pg,rec_Dt/alpha);
      mat.ir(2,k).ir(4,ndim+k).axpy(mass_pg,-1.0);
    }
    mat.rs();

  }
  shape.rs();
  res.rs();
  // tmp4.ctr(mat,2,1,4,3);
  // tmp4.print(nel*ndof);
    
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_load
::init() {
  
#define MAXPROPS 100
  elprpsindx.mono(MAXPROPS);
  propel.mono(MAXPROPS);
  
  // TGETOPTNDEF_ND(thash,int,ndim,none); //nd

  int ierr, iprop=0;
  pressure_indx = iprop; 
  ierr = get_prop(iprop,elem_prop_names,
		  thash,elprpsindx.buff(),propel.buff(), 
		  "pressure",1);
  nprops = iprop;
  nor.resize(1,ndim);
  Jaco.resize(2,ndimel,ndim);
  tmp.resize(2,nel,ndim);
}

#define ELEMPROPS(j,k) VEC2(elemprops,j,k,nelprops)

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void ld_elasticity_load
::element_connector(const FastMat2 &xloc,
		    const FastMat2 &state_old,
		    const FastMat2 &state_new,
		    FastMat2 &res,FastMat2 &mat){
  res.set(0.);
  mat.set(0.);
  load_props(propel.buff(),elprpsindx.buff(),nprops,
	     &(ELEMPROPS(elem,0)));
  double pressure = *(propel.buff()+pressure_indx);
  res.is(2,ndim+1,2*ndim);

  xstar.set(xloc);

  state.set(state_old);
  state.is(2,1,ndim);
  xstar.axpy(state,1-alpha);
  state.rs();

  state.set(state_new);
  state.is(2,1,ndim);
  xstar.axpy(state,alpha);
  state.rs();

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
    tmp.prod(shape,nor,1,2);
    res.axpy(tmp,-pressure);
#if 0
    printf("elem %d, wpgdet %f\n",elem,wpgdet);
    nor.print("nor: ");
#endif
  }
  dshapexi.rs();
  shape.rs();
  res.rs();
}

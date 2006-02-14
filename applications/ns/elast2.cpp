//__INSERT_LICENSE__
//$Id: elast2.cpp,v 1.8 2006/02/14 23:48:30 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "nsgath.h"
#include "elast2.h"

void elasticity2::init() {

  int ierr;
  //o Young modulus
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);

  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);

  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;
  // Dos opciones para imprimir
  // printf("rec_Dt: %d\n",rec_Dt);
  // SHV(rec_Dt);
  assert(!(rec_Dt>0. && rho==0.));
  assert(ndof==2*ndim);

  ntens = ndim*(ndim+1)/2;
  nen = nel*ndim;
  
  B.resize(2,ntens,nen).set(0.);
  C.resize(2,ntens,ntens).set(0.);
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  mass_pg.resize(2,nel,nel);

  // Plane strain
  if (ndim==2) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), c2=E/(2.*(1.+nu)),
      c3=nu/(1.-nu);
    C.setel(c1,1,1)
      .setel(c1*c3,1,2)
      .setel(c1*c3,2,1)
      .setel(c1,2,2)
      .setel(c2,3,3);
#if 0
    // From RAPPORT SF-101, LTAS, Godinas A., Jetteur P., Laschet G.
    // This modification of the constitutive relation should avoid
    // blockage of thin elements for plate model. 
    int plate_model = 1;
    if (plate_model) {
      printf("using plate model\n");
      C.set(0.0);
      double c1 = E/(1.0-nu*nu);
      C.setel(c1,1,1);
      C.setel(E,2,2);
      C.setel(c1*(5.0/6.0)*(1.0-nu)/2.0,3,3);
    }
#endif
  } else if (ndim==3) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), 
      c2 = (1-2*nu)/2./(1-nu),
      c3=nu/(1.-nu);
      C.is(1,1,3).is(2,1,3).set(c3)
	.rs().d(1,2).is(1,1,3).set(1.)
	.rs().d(1,2).is(1,4,6).set(c2)
	.rs().scale(c1);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"wrong dimension: %d\n",ndim);
    assert(0);
  }
}

void elasticity2::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  B.reshape(3,ntens,nel,ndim);
  res.set(0.0);
  mat.set(0.0);

  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
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
    
    // construccion de matriz B
    if (ndim==2) {
      B.ir(1,1).ir(3,1).set(dshapex.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex.ir(1,2));
      B.ir(1,3).ir(3,1).set(dshapex.ir(1,2));
      B.ir(1,3).ir(3,2).set(dshapex.ir(1,1));
    } else if (ndim==3) {
      B.ir(1,1).ir(3,1).set(dshapex.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex.ir(1,2));
      B.ir(1,3).ir(3,3).set(dshapex.ir(1,3));

      B.ir(1,4).ir(3,2).set(dshapex.ir(1,1));
      B.ir(1,4).ir(3,1).set(dshapex.ir(1,2));

      B.ir(1,5).ir(3,3).set(dshapex.ir(1,1));
      B.ir(1,5).ir(3,1).set(dshapex.ir(1,3));

      B.ir(1,6).ir(3,3).set(dshapex.ir(1,2));
      B.ir(1,6).ir(3,2).set(dshapex.ir(1,3));
    }
    B.rs();
    dshapex.rs();
    res.rs();

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

    xstar.set(xnew).scale(alpha).axpy(xold,1-alpha);
    vstar.set(vnew).scale(alpha).axpy(vold,1-alpha);

    strain.prod(B,xstar,1,-1,-2,-1,-2);
    stress.prod(C,strain,1,-1,-1);

    shape.ir(2,ipg+1);

    // Eq. for velocity (dofs=[ndim+1,2*ndim])
    // Res = [ -M*(xnew-xold)/dt + M*vstar;
    //         -M*rho*(vnew-vold)/dt - K*xstar
    // Jac = -D(res)/DU = [M/(Dt*alpha) -M;
    //                     K             M*rho/(Dt*alpha)]
    // (recordar que Jac debe ser -1/alpha*d(Res)/d(X) 

    // Inertia term
    a.set(vnew).rest(vold);
    tmp.prod(shape,a,-1,-1,1);
    tmp2.prod(shape,tmp,1,2);
    res.is(2,ndim+1,2*ndim).axpy(tmp2,-wpgdet*rec_Dt*rho);

    // Elastic force residual computation
    res_pg.prod(B,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet).rs();
    
    mass_pg.prod(shape,shape,1,2).scale(wpgdet*rec_Dt*rho/alpha);
    for (int k=1; k<=ndim; k++) 
      mat.ir(2,ndim+k).ir(4,ndim+k).add(mass_pg);
    mat.rs();

    // Jacobian computation
    mat_pg1.prod(C,B,1,-1,-1,2,3);
    mat_pg2.prod(B,mat_pg1,-1,1,2,-1,3,4);
    mat.is(2,ndim+1,2*ndim).is(4,1,ndim)
      .axpy(mat_pg2,wpgdet).rs();
    
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
void elast_energy_integrator::init(){
  int ierr;
  assert(gather_length==1 || gather_length==2);
  TGETOPTDEF(thash,double,Young_modulus,0.);
  E=Young_modulus;
  assert(Young_modulus>0.);
  //o Poisson ratio
  TGETOPTDEF(thash,double,Poisson_ratio,0.);
  nu=Poisson_ratio;
  assert(nu>=0. && nu<0.5);
  //o Density
  TGETOPTDEF(thash,double,density,0.);
  rho=density;
  assert(rho>0.);
  //o ndim
  TGETOPTDEF_ND(thash,int,ndim,0);
  assert(ndim>0);
  assert(ndof==2*ndim);

  TGETOPTDEF(thash,int,ndimel,ndim);
  
  e_total.resize(0).set(0.);
  e_kin.resize(0).set(0.);
  e_pot.resize(0).set(0.);
  ntens = ndim*(ndim+1)/2;
  C.resize(2,ntens,ntens).set(0.);
  stress.resize(1,ntens);
  strain.resize(1,ntens);
  // Plane strain
  if (ndim==2) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), c2=E/(2.*(1.+nu)),
      c3=nu/(1.-nu);
    C.setel(c1,1,1)
      .setel(c1*c3,1,2)
      .setel(c1*c3,2,1)
      .setel(c1,2,2)
      .setel(c2,3,3);
  } else if (ndim==3) {
    double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), 
      c2 = (1-2*nu)/2./(1-nu),
      c3=nu/(1.-nu);
    C.is(1,1,3).is(2,1,3).set(c3)
      .rs().d(1,2).is(1,1,3).set(1.)
      .rs().d(1,2).is(1,4,6).set(c2)
      .rs().scale(c1);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"wrong dimension: %d\n",ndim);
    assert(0);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void elast_energy_integrator
::set_pg_values(vector<double> &pg_values,FastMat2 &u,
		FastMat2 &uold,FastMat2 &xpg,FastMat2 &Jaco,
		double wpgdet,double time) {
  int ierr;
  if (ndim==2){
    double eps_xx = grad_u.get(1,1);
    strain.setel(eps_xx,1);
    double eps_yy = grad_u.get(2,2);
    strain.setel(eps_yy,2);
    double eps_xy = (grad_u.get(1,2) + grad_u.get(2,1));
    strain.setel(eps_xy,3);
    stress.prod(C,strain,1,-1,-1);
    e_pot.prod(strain,stress,-1,-1).scale(0.5);
  } else if (ndim==3){
    double eps_xx = grad_u.get(1,1); strain.setel(eps_xx,1);
    double eps_yy = grad_u.get(2,2); strain.setel(eps_yy,2);
    double eps_zz = grad_u.get(3,3); strain.setel(eps_zz,3);
    double eps_xy = (grad_u.get(1,2) + grad_u.get(2,1));
    double eps_xz = (grad_u.get(1,3) + grad_u.get(3,1));
    double eps_yz = (grad_u.get(2,3) + grad_u.get(3,2));
    strain.setel(eps_xy,4);
    strain.setel(eps_xz,5);
    strain.setel(eps_yz,6);
    stress.prod(C,strain,1,-1,-1);
    e_pot.prod(strain,stress,-1,-1).scale(0.5);
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"wrong dimension: %d\n",ndim);
    assert(0);
  }
  tmp1.set(u.is(1,ndim+1,ndof));
  u.rs();
  e_kin.prod(tmp1,tmp1,-1,-1).scale(rho).scale(0.5);
  e_total.set(e_pot).add(e_kin);  
  
  if (gather_length==1) {
    pg_values[0] = double(e_total)*wpgdet;    
  } else {
    pg_values[0] = double(e_kin)*wpgdet;    
    pg_values[1] = double(e_pot)*wpgdet;    
  } 
}

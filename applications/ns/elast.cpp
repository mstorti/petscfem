//__INSERT_LICENSE__
//$Id: elast.cpp,v 1.13 2003/09/17 00:51:31 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include "./fm2funm.h"

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

extern int SIZE, MY_RANK;

class MyFun2 : public FastMat2_funm {
public:
  // Some kind of `hardening'
  // double f(double l) { return pow(l,-0.3); }

  // No hardening at all
  double f(double l) { return 1.; }
} my_fun2;

#if 0
class MyFun : public FastMat2_fund {
public:
  void f(const FastMat2 &L,FastMat2 &fL) { fL.set(L); }
} my_fun;
#endif

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
class my_locker {
public:
  int cookie;
  double dcookie;
};

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#define M 20
void elasticity::element_init() {
  my_locker *locker = new my_locker;
  if (!locker) {
    set_error(1);
    return;
  }
  locker->cookie = rand();
  locker->dcookie = drand();
  local_store_address(elem) = locker;
  if (!(elem % M)) printf("[%d] element %d, cookie %d, dcookie %g\n",
			  MY_RANK, elem, locker->cookie, locker->dcookie);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void elasticity::init() {

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
  assert(ndof==ndim);

  ntens = ndim*(ndim+1)/2;
  nen=nel*ndof;
  
  // tal vez el resize blanquea
  B.resize(2,ntens,nen).set(0.);
  C.resize(2,ntens,ntens).set(0.);
  Jaco.resize(2,ndim,ndim);
  G.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  

  my_fun2.init(G);
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

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  my_locker *locker = (my_locker *)local_store_address(elem);
  if (!(elem % M)) printf("[%d] At element_connector:start: "
			  "element %d, cookie %d, dcookie %g\n",
			  MY_RANK, elem, locker->cookie, locker->dcookie);
  locker->cookie = int(10000*rand());
  locker->dcookie = drand();
  if (!(elem % M)) printf("[%d] In element_connector, setting: "
			  "element %d, cookie %d, dcookie %g\n",
			  MY_RANK, elem, locker->cookie, locker->dcookie);

  B.reshape(3,ntens,nel,ndof);

  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
    // x_new.set(xloc).add(state_new);
    x_new.set(xloc);
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    // Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    Jaco.prod(dshapexi,x_new,1,-1,-1,2);
    G.prod(Jaco,Jaco,-1,1,-1,2);
    my_fun2.apply(G,fG);
    
    double detJaco = Jaco.det();
    if (detJaco<=0.) {
      detj_error(detJaco,elem);
      set_error(1);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);
    iJaco.inv(Jaco);
    dshapex.prod(iJaco,dshapexi,1,-1,-1,2);
    dshapex_scaled.prod(fG,dshapex,1,-1,-1,2);
    //dshapex_scaled.set(dshapex);
    
    // Recall: \epsilon = B dudx
    // where \epsilon = [e_xx e_yy e_xy] (2D)
    //                  [e_xx e_yy e_zz e_xy e_xz e_yz] (3D)
    // dudx is de gradient of displacements. 
    if (ndim==2) {
      B.ir(1,1).ir(3,1).set(dshapex_scaled.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,1).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,2).set(dshapex_scaled.ir(1,1));
    } else if (ndim==3) {
      B.ir(1,1).ir(3,1).set(dshapex_scaled.ir(1,1));
      B.ir(1,2).ir(3,2).set(dshapex_scaled.ir(1,2));
      B.ir(1,3).ir(3,3).set(dshapex_scaled.ir(1,3));

      B.ir(1,4).ir(3,2).set(dshapex_scaled.ir(1,1));
      B.ir(1,4).ir(3,1).set(dshapex_scaled.ir(1,2));

      B.ir(1,5).ir(3,3).set(dshapex_scaled.ir(1,1));
      B.ir(1,5).ir(3,1).set(dshapex_scaled.ir(1,3));

      B.ir(1,6).ir(3,3).set(dshapex_scaled.ir(1,2));
      B.ir(1,6).ir(3,2).set(dshapex_scaled.ir(1,3));
    } else assert(0);

    dshapex_scaled.rs();
    
    // B.rs().reshape(2,ntens,nen);
    B.rs();
    strain.prod(B,state_new,1,-1,-2,-1,-2);
    stress.prod(C,strain,1,-1,-1);
    
    // Residual computation
    res_pg.prod(B,stress,-1,1,2,-1);
    res.axpy(res_pg,-wpgdet);
    
    // Jacobian computation
    mat_pg1.prod(C,B,1,-1,-1,2,3);
    mat_pg2.prod(B,mat_pg1,-1,1,2,-1,3,4);
    mat.axpy(mat_pg2,wpgdet);
    
  }
    
}

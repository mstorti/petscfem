//__INSERT_LICENSE__
//$Id: elast.cpp,v 1.6 2002/12/09 02:57:53 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

class FastMat2_fund {
private:
  int m;
  FastMat2 V,D,tmp,lambda,flambda;
public:
  void init(FastMat2 &A);
  void apply(const FastMat2 &A,FastMat2 &fA);
  virtual void f(const FastMat2 &D,FastMat2 &fD)=0;
};

void FastMat2_fund::init(FastMat2 &A) {
  assert(A.n()==2);
  assert(A.dim(1)==A.dim(2));
  m = A.dim(1);
  V.resize(2,m,m);
  D.resize(2,m,m).set(0.);
  lambda.resize(1,m);
  flambda.resize(1,m);
  tmp.resize(2,m,m);
}

void FastMat2_fund::apply(const FastMat2 &A,FastMat2 &fA) {
  lambda.seig(A,V);
  f(lambda,flambda);
  D.d(1,2).set(flambda).rs();
  tmp.prod(D,V,1,-1,2,-1);
  fA.prod(V,tmp,1,-1,-1,2);
}

double FastMat2_funm_ff(double x,void *a);

class FastMat2_funm : public FastMat2_fund {
private:
  void f(const FastMat2 &D,FastMat2 &fD);
public:
  virtual double f(double)=0;
};

double FastMat2_funm_ff(double x,void *a) {
  FastMat2_funm *fff = (FastMat2_funm *)a;
  return fff->f(x);
}

void FastMat2_funm::f(const FastMat2 &D,FastMat2 &fD) {
  fD.set(D).fun(FastMat2_funm_ff,this); 
}

class MyFun2 : public FastMat2_funm {
  // double f(double l) { 1./sqrt(l); }
  double f(double l) { return l; }
} my_fun2;

class MyFun : public FastMat2_fund {
  void f(const FastMat2 &L,FastMat2 &fL) { fL.set(L); }
} my_fun;

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

void elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  B.reshape(3,ntens,nel,ndof);

  // loop over Gauss points
  for (int ipg=0; ipg<npg; ipg++) {
    
    x_new.set(xloc).add(state_new);
    dshapexi.ir(3,ipg+1); // restriccion del indice 3 a ipg+1
    // Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    Jaco.prod(dshapexi,x_new,1,-1,-1,2);
    G.prod(Jaco,Jaco,-1,1,-1,2);
    my_fun2.apply(G,fG);
    
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
    dshapex_scaled.prod(fG,dshapex,1,-1,-1,2);
    
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

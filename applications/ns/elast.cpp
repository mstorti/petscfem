//__INSERT_LICENSE__
//$Id: elast.cpp,v 1.2 2001/11/30 12:41:08 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

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

  // 2D
  ntens=3;
  nen=nel*ndof;
  
  // tal vez el resize blanquea
  B.resize(2,ntens,nen).set(0.);
  C.resize(2,ntens,ntens).set(0.);
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  

  // Plane strain
  double c1=E*(1.-nu)/((1.+nu)*(1.-2.*nu)), c2=E/(2.*(1.+nu)),
         c3=nu/(1.-nu);
  C.setel(c1,1,1).setel(c1*c3,1,2)
    .setel(c1*c3,2,1).setel(c1,2,2)
    .setel(c2,3,3);

}

void elasticity::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  B.reshape(3,ntens,nel,ndof);

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
    B.ir(1,1).ir(3,1).set(dshapex.ir(1,1));
    B.ir(1,2).ir(3,2).set(dshapex.ir(1,2));
    B.ir(1,3).ir(3,1).set(dshapex.ir(1,2));
    B.ir(1,3).ir(3,2).set(dshapex.ir(1,1));

    dshapex.rs();
    
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

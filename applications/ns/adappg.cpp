//__INSERT_LICENSE__
//$Id: adappg.cpp,v 1.1 2001/11/30 13:22:24 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "elast.h"

void adaptor_pg::init() {
  nen=nel*ndof;
  Jaco.resize(2,ndim,ndim);
  dshapex.resize(2,ndim,nel);  
  grad_state_new.resize(2,ndim,ndof);
  grad_state_old.resize(2,ndim,ndof);

}

void adaptor_pg::clean() {
  nen=nel*ndof;
  Jaco.resize();
  dshapex.resize();  
  grad_state_new.resize();
  grad_state_old.resize();

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

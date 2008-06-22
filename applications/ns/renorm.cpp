//__INSERT_LICENSE__
//$Id: elast2.cpp,v 1.11.20.1 2007/02/19 20:23:56 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "renorm.h"

void renorm::init() {
  int ierr;

  PETSCFEM_ASSERT0(ndof==1,"Only ndof==1 considered");  
  PETSCFEM_ASSERT0(ndim==1,"Only ndim==1 considered");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,creac,NAN);
  PETSCFEM_ASSERT0(!isnan(creac),"creac is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,kond,NAN);
  PETSCFEM_ASSERT0(!isnan(kond),"kond is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,mpenal,NAN);
  PETSCFEM_ASSERT0(!isnan(mpenal),"mpenal is required");  

}

void renorm::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  res.set(0.0);
  mat.set(0.0);

  // loop over Gauss points
  phi.set(state_new);
  phi.ir(2,1);
  res.ir(2,1);
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
    
    grad_phi.prod(phi,dshapex,-1,1,-1);
    tmp.prod(dshapex,grad_phi,-1,1,-1);
    res.axpy(tmp,kond*wpgdet);

    shape.ir(2,ipg+1);
    phipgm.prod(phi,shape,-1,-1);
    double phipg = phipgm; 
    tmp2.set(shape).scale(phipg);
    res.axpy(tmp2,creac*wpgdet);
  }
  phi.rs();
  shape.rs();
  res.rs();
  // tmp4.ctr(mat,2,1,4,3);
  // tmp4.print(nel*ndof);
    
}


//__INSERT_LICENSE__
//$Id: elast2.cpp,v 1.11.20.1 2007/02/19 20:23:56 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "renorm2.h"

void renorm2::init() {
  int ierr;

  PETSCFEM_ASSERT0(ndof==1,"Only ndof==1 is considered");  

  //o Constant affecting #||\grad phi||^2# term
  TGETOPTDEF_ND(thash,double,c_grad_phi,NAN);
  PETSCFEM_ASSERT0(!isnan(c_grad_phi),"c_grad_phi is required");  

  //o Constant affecting Laplacian term
  TGETOPTDEF_ND(thash,double,kond,NAN);
  PETSCFEM_ASSERT0(!isnan(kond),"kond is required");  

  //o Constant affecting reaction term
  TGETOPTDEF_ND(thash,double,c_reac,NAN);
  PETSCFEM_ASSERT0(!isnan(c_reac),
                   "c_reac is required");  

  //o Constant affecting reaction term
  TGETOPTDEF_ND(thash,int,use_lumped,0);
}

void renorm2::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  res.set(0.0);
  mat.set(0.0);

  // loop over Gauss points
  phi.set(state_new);
  phi.ir(2,1);
  res.ir(2,1);
  double vol=0.0;
  for (int ipg=0; ipg<npg; ipg++) {
    dshapexi.ir(3,ipg+1);
    shape.ir(2,ipg+1);
    Jaco.prod(dshapexi,xloc,1,-1,-1,2);
    
    double detJaco = Jaco.det();
    if (detJaco <= 0.) {
      printf("Jacobian of element %d is negative or null\n"
	     " Jacobian: %f\n",elem,detJaco);
      PetscFinalize();
      exit(0);
    }
    double wpgdet = detJaco*wpg.get(ipg+1);
    vol += wpgdet;
    iJaco.inv(Jaco);
    dshapex.prod(iJaco,dshapexi,1,-1,-1,2);
    
    phipgm.prod(phi,shape,-1,-1);
    double phipg = phipgm; 
    double fphi = phipg*(phipg*phipg-1.0);

    grad_phi.prod(phi,dshapex,-1,1,-1);
    double gphi2 = grad_phi.sum_square_all();
    // double q = fphi*c_reac - c_grad_phi*phipg*sqrt(gphi2);
    double q = (!use_lumped)*fphi*c_reac 
      + c_grad_phi*phipg*sqrt(gphi2);
    res.axpy(shape,q*wpgdet);

    tmp2.prod(dshapex,grad_phi,-1,1,-1);
    res.axpy(tmp2,kond*wpgdet);
  }
  if (use_lumped) {
    double 
      *phip = phi.storage_begin(),
      *resp = res.storage_begin(),
      volnod = vol/nel;
    // printf("volnod %f\n",volnod);
    for (int j=0; j<nel; j++) {
      double phix = phip[j];
      double fphi = phix*(phix*phix-1.0);
      // printf("phi %f, fphi %f\n",phix,fphi);
      resp[j] += volnod*fphi;
    }
  }

  phi.rs();
  shape.rs();
  res.rs();
}

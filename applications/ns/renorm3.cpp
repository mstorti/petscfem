//__INSERT_LICENSE__
//$Id: elast2.cpp,v 1.11.20.1 2007/02/19 20:23:56 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "renorm3.h"

extern double total_liquid_volume;

void renorm3::init() {
  int ierr;

  PETSCFEM_ASSERT0(ndof==1,"Only ndof==1 considered");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,creac,NAN);
  PETSCFEM_ASSERT0(!ISNAN(creac),"creac is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,kond,NAN);
  PETSCFEM_ASSERT0(!ISNAN(kond),"kond is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,mpenal,NAN);

  //o Constant affecting reaction term
  TGETOPTDEF_ND(thash,int,use_lumped,0);

}

void renorm3::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  res.set(0.0);
  mat.set(0.0);

  // loop over Gauss points
  phi.set(state_new);
  phi.ir(2,1);
  phiold.set(state_old);
  phiold.ir(2,1);
  res.ir(2,1);
  vel = 0; // elemental volume
  double vol=0.0;
  for (int ipg=0; ipg<npg; ipg++) {
    dshapexi.ir(3,ipg+1);
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
    double fphi = phipg*(phipg*phipg-1.0);
    tmp2.set(shape).scale(fphi);
    res.axpy(tmp2,creac*wpgdet*(!use_lumped));

    // penality term by continuous function
    phipgmold.prod(phiold,shape,-1,-1);
    double phipgold = phipgmold; 
    // tanh 
    double pi = 4*atan(1.0);
    double tphi = tanh(2*pi*phipg);
    tmp3.set(shape).scale(tphi);
    res.axpy(tmp3,mpenal*wpgdet);
    double tphiold = tanh(2*pi*phipgold);
    tmp4.set(shape).scale(tphiold);
    res.axpy(tmp4,-mpenal*wpgdet);

    // volume(area) of the phi>0 region
    //    int liq = 1;
    int liq = phipg>0;
    vol += wpgdet;
    vel += wpgdet*liq;
  }
  total_liquid_volume += vel/(npg+1);
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
      resp[j] += creac*volnod*fphi;
    }
  }
  
   phi.rs();
   phiold.rs();
   shape.rs();
   res.rs();
   xlocc.rs();
    
}


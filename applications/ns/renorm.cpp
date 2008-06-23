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
  PETSCFEM_ASSERT0(ndim==1 || ndim==2,
                   "Only ndim==1,2 implemented so far");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,creac,NAN);
  PETSCFEM_ASSERT0(!isnan(creac),"creac is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,kond,NAN);
  PETSCFEM_ASSERT0(!isnan(kond),"kond is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,mpenal,NAN);
  PETSCFEM_ASSERT0(!isnan(mpenal) && mpenal>0 
                   && nel==2 && ndim==1,
                   "mpenal only implemented for linear segements");  

  resh.resize(1,nel);
  C.resize(2,nel,nel);
  phirot.resize(1,nel);
  xrot.resize(2,nel,ndim);
}

void renorm::compute_H_term(FastMat2 &phi) {
  double phimax = phi.max_all();
  double phimin = phi.min_all();
  resh.set(0.0);
#if 1
  if ((phimax>0.0) != (phimin>0.0)) {
    double 
      *phip = phi.storage_begin(),
      *reshp = resh.storage_begin();
    if (ndim==1) {
      double 
        x1 = xlocc.get(1,1),
        x2 = xlocc.get(2,1),
        h = fabs(x2-x1), hliq=NAN, xiav=1.0;
      double
        phi1 = phip[0],
        phi2 = phip[1],
        xii = -phi1/(phi2-phi1),
        xi1, xi2;
      if (phi1>0) { xi1 = 0; xi2 = xii; }
      else { xi1 = xii; xi2 = 1.0; }
        
      xiav = (xi1+xi2)/2;
      hliq = h*(xi2-xi1);
      reshp[0] = (1-xiav)*hliq;
      reshp[1] = xiav*hliq;
    } else if (ndim==2) {
      // Nbr of positive values
      int pos = 0;
      for (int j=0; j<nel; j++)
        pos += phip[j]>0.0;
      assert(pos==1 || pos==2);
      // If 2 positive and 1 negative we reduce
      // to the 1 positive case and after rest the
      // entalphy of the whole element
      // flag := indicates wether we performed the
      // inversion or not
      flag = 0;
      if (pos==2) {
        for (int j=0; j<nel; j++) 
          phip[j] = -phip[j];
        flag = 1;
      }
      // Indx of vertex with positive phi
      int pvrtx=-1;
      for (int j=0; j<nel; j++) {
        if (phip[j]>0.0) pvrtx=j;
        break;
      }
      assert(pvrtx>=0);
      // Rotate phi and coords
      double
        *phirotp = phirot.storage_begin(),
        *xrotp = xrot.storage_begin(),
        *xlocp = xlocc.storage_begin();
      for (int j=0; j<nel; j++) {
        int jj = modulo(j-pvrtx,nel);
        phirotp[jj] = phip[j];
        for (int k=0; k<ndim; k++) 
          xrotp[j*ndim+k] = xlocp[j*ndim+k];
      }
    }
  } 
  // else hliq = h*(phimax>0.0);
#endif
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
  phiold.set(state_old);
  phiold.ir(2,1);
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
    double fphi = phipg*(phipg*phipg-1.0);
    tmp2.set(shape).scale(fphi);
    res.axpy(tmp2,creac*wpgdet);

    if (!isnan(mpenal) && mpenal>0) {
      xlocc.set(xloc);
      compute_H_term(phi);
      res.axpy(resh,mpenal);
      compute_H_term(phiold);
      res.axpy(resh,-mpenal);
    }
  }
  phi.rs();
  phiold.rs();
  shape.rs();
  res.rs();
  // tmp4.ctr(mat,2,1,4,3);
  // tmp4.print(nel*ndof);
    
}


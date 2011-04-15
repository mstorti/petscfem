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

extern double total_liquid_volume;

void renorm::init() {
  int ierr;

  PETSCFEM_ASSERT0(ndof==1,"Only ndof==1 considered");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,creac,NAN);
  PETSCFEM_ASSERT0(!isnan(creac),"creac is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,kond,NAN);
  PETSCFEM_ASSERT0(!isnan(kond),"kond is required");  

  //o Poisson ratio
  TGETOPTDEF_ND(thash,double,mpenal,NAN);

  resh.resize(1,nel);
  C.resize(2,nel,nel);
  phirot.resize(1,nel);
  xrot.resize(2,nel,ndim);
}

void renorm::compute_H_term(FastMat2 &phi) {
  double phimax = phi.max_all();
  double phimin = phi.min_all();
  resh.set(0.0);
  double 
    *phip = phi.storage_begin(),
    *reshp = resh.storage_begin();
  if (ndim==1) {
    double 
      x1 = xlocc.get(1,1),
      x2 = xlocc.get(2,1),
      h = fabs(x2-x1);
    if ((phimax>0.0) != (phimin>0.0)) {
      assert(nel==2);
      double hliq=NAN, xiav=1.0;
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
    } else if (phimin>0.0) {
      reshp[0] = h/nel;
      reshp[1] = h/nel;
    }
  } else if (ndim==2) {
    assert(nel==3);
    if ((phimax>0.0) != (phimin>0.0)) {
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
        if (phip[j]>=0.0) {
	  pvrtx=j;
	  break;
	}
      }
      assert(pvrtx>=0);
      // Rotate phi and coords
      double
        *phirotp = phirot.storage_begin(),
        *xrotp = xrot.storage_begin(),
        *xlocp = xlocc.storage_begin(),
        *xareap = xarea.storage_begin();
      int idxrot[nel];
      //      for (int j=0;j<nel;j++) idxrot[j] = j;
      for (int j=0; j<nel; j++) {
        int jj = modulo(j-pvrtx,nel);
        phirotp[jj] = phip[j];
	idxrot[j]=jj;
        for (int k=0; k<ndim; k++) 
          xrotp[jj*ndim+k] = xlocp[j*ndim+k];
        for (int k=0; k<ndim; k++) 
          xareap[jj*ndim+k] = xlocp[j*ndim+k];
      }
      // Here we have in xrotp and phirotp
      // a state element with the positive value
      // in the first node 
      /* --------- */
      // Auxiliar points for determining the 
      // triangular area 
      // First, location of points phi=0 over local edges
      // auxiliar nodal coordinates and center of gravity 
      //      xarea.set(xrot);
      double prop[nel-1];
      for (int k=0; k<nel-1; k++) {
	double deno = phirotp[k+1]-phirotp[0];
	if (abs(deno)>0.000001) {
	  prop[k] = -phirotp[0]/deno;
	} else {
	  prop[k] = 1;
	}
       	if (abs(prop[k])>1.0) 
	  printf("prop[k] %f phirotp[0] %f phirotp[k+1] %f\n",prop[k],phirotp[0],phirotp[k+1]);
	assert(prop[k]<=1.0 && prop[k]>=0.0);
	for (int j=0; j<ndim; j++){
	  xareap[ndim*(k+1)+j] = xrotp[j]+prop[k]*(xrotp[ndim*(k+1)+j]-xrotp[j]);
	}
       }
      //      xarea.print();
      // Triangular area
      t_area=0.0;
      t_area = fabs(0.5*(xareap[0]*xareap[3]-xareap[3]*xareap[4]-xareap[0]*xareap[5]
			 -xareap[1]*xareap[2]+xareap[5]*xareap[2]+xareap[4]*xareap[1]));
      // Nodal resh
      // Values of shape functions over nodes of the subtriangle 
      double reshrot[nel];
      double scale = t_area/3;
      reshrot[0] = scale*(3-prop[0]-prop[1]);
      reshrot[1] = scale*(prop[0]);
      reshrot[2] = scale*(prop[1]);
      // Procedure for subtracting area H(phi)=0;
      inv = 1; 
      for (int j=0; j<nel; j++) reshp[j] = 0.0;
      double area=NAN;
      if (flag) {
	inv = -1;
	area = fabs(0.5*(xrotp[0]*xrotp[3]-xrotp[3]*xrotp[4]-xrotp[0]*xrotp[5]
			 -xrotp[1]*xrotp[2]+xrotp[5]*xrotp[2]+xrotp[4]*xrotp[1]));
	for (int j=0; j<nel; j++) reshp[j] = area/3;
      }
      // Values resh from reshrot
      for (int j=0; j<nel; j++){
	int jj = idxrot[j];
	reshp[j] += inv*reshrot[jj];      
	if (reshp[j]<0) printf("reshp[j] %f reshrot[j] %f area %f t_area %f reshsum %f\n",
			       reshp[j],reshrot[jj],area,t_area,reshrot[0]+reshrot[1]+reshrot[2]);
	assert(reshp[j]>=0.0);
      }
    } else if (phimin>0.0) {
      double *xlocp = xlocc.storage_begin();
      double area = fabs(0.5*(xlocp[0]*xlocp[3]-xlocp[3]*xlocp[4]-xlocp[0]*xlocp[5]
			 -xlocp[1]*xlocp[2]+xlocp[5]*xlocp[2]+xlocp[4]*xlocp[1]));
      for (int j=0; j<nel; j++) reshp[j] = area/3;
    } else if (phimax<0.0) {//creo que no hace falta
      for (int j=0; j<nel; j++) reshp[j] = 0.0;
    }
  }
  phirot.rs();
}

void renorm::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  res.set(0.0);
  resh.set(0.0);
  mat.set(0.0);

  // loop over Gauss points
  phi.set(state_new);
  phi.ir(2,1);
  phiold.set(state_old);
  phiold.ir(2,1);
  res.ir(2,1);
  vel = 0; // elemental volume
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
    // inicio correccion
    //    if (phipg>1.10) phipg=1.10;
    //if (phipg<-1.10) phipg=-1.10;
    // fin correccion
    double fphi = phipg*(phipg*phipg-1.0);
    tmp2.set(shape).scale(fphi);
    res.axpy(tmp2,creac*wpgdet);

    // penality term by continuous function
    // mpenal*(phi-phiold)
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
    vel += wpgdet*liq;
    //    printf("liq %d \n",liq);
    //    printf("vel %f \n",vel);
  }
 total_liquid_volume += vel/(npg+1);

  //  printf("total mass %f \n",total_liquid_volume);

//   if (!isnan(mpenal) && mpenal>0) {
//     xlocc.set(xloc);
//     xarea.set(xloc);
//     compute_H_term(phi);
//     res.axpy(resh,mpenal);
//     compute_H_term(phiold);
//     res.axpy(resh,-mpenal);
//   }
//   //  resh.print();
//   total_vol = resh.sum_all();

  phi.rs();
  phiold.rs();
  shape.rs();
  res.rs();
  xlocc.rs();
  xarea.rs();
  resh.rs();
  // tmp4.ctr(mat,2,1,4,3);
  // tmp4.print(nel*ndof);
    
}


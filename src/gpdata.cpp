/*
  This file belongs to the PETSc - FEM package, a library and
  application suite oriented to the Finite Element Method based on PETSc. 
  Copyright (C) 1999, 2000  Mario Alberto Storti
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.

*/

#include "sles.h"
#include <math.h>
#include <fnmatch.h>

#include "fem.h"
#include "util2.h"
#include "gpdata.h"

#define GPERROR \
    {PFEM_TRACE(""); \
    PetscPrintf(PETSC_COMM_WORLD,"Not implemented combination: geometry=\"%s\","\
           "nel=%d, npg=%d, ndimel=%d, \n", \
	   geom,nel,npg,ndimel); \
    PetscFinalize(); \
		       exit(0);} \

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
#undef __FUNC__
#define __FUNC__ "cartesian_2d_shape(Rowvector &,Matrix &,matrix &)"
int cartesian_2d_shape(RowVector &shape,Matrix &dshapexi,
		       double xipg,double etapg) {

  int nel=4,ndim=2;
  static int flag=0;
  static Matrix xinode;
  if (flag==0) {
    flag=1;
    xinode.ReSize(2,4);
    xinode << -1 << 1 << 1 << -1 << -1 << -1 << 1 << 1;
  }

  shape.ReSize(nel);
  dshapexi.ReSize(ndim,nel);

  for (int iloc=1; iloc<=nel; iloc++) {
    // 1D Shape Functions
    // along X
    double sxi=(1.+xipg*xinode(1,iloc))/2.;
    double dsxidxi=xinode(1,iloc)/2.;

    // along Y
    double seta=(1.+etapg*xinode(2,iloc))/2.;
    double dsetadeta=xinode(2,iloc)/2.;
	  
    shape(iloc)=sxi*seta;
    dshapexi(1,iloc)=dsxidxi*seta;
    dshapexi(2,iloc)=sxi*dsetadeta;
  }
  return 0;
}


GPdata::GPdata(const char *geom,int ndimel,int nel,int npg_,int
	       mat_version_=GP_NEWMAT) {
  mat_version = mat_version_;
  npg= npg_;
  int ipg;
  // Dimension GP vectors

  // ndimel:= dimension of the element. Cartesian elements are the
  //          cartesian product of 1D linear elements. ndimel is not
  //          necessarily equal to the number of dimensions in the
  //          space. For instance, a boundary element (for imposed
  //          flux or convection boundary condition) has a dimension
  //          lower than the spatial dimension.
  
  // Shape function
  shape = new RowVector[npg];
  // Weights
  wpg = new double[npg];
  // dshapexi:= [ipg](jd,jloc) Gradient of shape function
  //    with respect to master element coordinates
  //    ipg = number of Gauss Point
  //    jd = spatial coordinate
  //    jloc = local node number
  dshapexi = new Matrix[npg];
  // dshapexi:= [ipg](jd,jloc) Gradient of shape function
  //    with respect to global coordinates
  dshapex = new Matrix[npg];

  if ( !(strcmp(geom,"triangle")) ) {

    double xipg,etapg;
    for (ipg=0; ipg<npg; ipg++) {
      if (npg==7) {
	// Three points in the edge centers
#define A1 (0.0597158717)
#define B1 (0.4701420641)
#define A2 (0.7974269853)
#define B2 (0.1012865073)

#define W1 (0.225)
#define W2 (0.1323941527)
#define W3 (0.1259391805)

	double wpg_[7]={W1,W2,W2,W2,W3,W3,W3};
	double xpg_[7][2]={
	  1./3., 1./3.,
	  A1, B1,
	  B1, A1,
	  B1, B1,
	  A2, B2,
	  B2, A2,
	  B2, B2};
	xipg  = xpg_[ipg][0];
	etapg = xpg_[ipg][1];
	wpg[ipg] = wpg_[ipg];

      } else if (npg==4) {
	// Three points in the edge centers
	double wpg_[4]={-27./48.,25./48.,25./48.,25./48.};
	double xpg_[4][2]={
	  1./3., 1./3.,
	  0.2, 0.2,
	  0.6, 0.2,
	  0.2, 0.6};
	xipg  = xpg_[ipg][0];
	etapg = xpg_[ipg][1];
	wpg[ipg] = wpg_[ipg];

      } else if (npg==3) {

	// Three points in the edge centers
	double a=0.5;
	if (ipg==0) {
	  xipg=1-2*a; etapg=a;
	} else if (ipg==1) {
	  xipg=a; etapg=1-2*a;
	} else {
	  xipg=a; etapg=a;
	}
	wpg[ipg] =  1/6.;

      } else if (npg==1) {
	// One point in the center of the element
	xipg=1./3.;
	etapg=1./3.;
	wpg[ipg] =  0.5;

      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Not valid value of npg. for triangles %d\n",
		    npg);
	CHKERRA(1); 
      }

      shape[ipg] = RowVector(nel);
      shape[ipg](1)=xipg;
      shape[ipg](2)=etapg;
      shape[ipg](3)=1-xipg-etapg;
    
      dshapex[ipg]= Matrix(ndimel,nel);
      dshapexi[ipg]= Matrix(ndimel,nel);
    
      dshapexi[ipg](1,1)=1.;
      dshapexi[ipg](2,1)=0.;
      dshapexi[ipg](1,2)=0.;
      dshapexi[ipg](2,2)=1.;
      dshapexi[ipg](1,3)=-1.;
      dshapexi[ipg](2,3)=-1.;
    }

  } else if ( !(strcmp(geom,"tetra")) ) {

    double xipg,etapg,zetapg;
    for (ipg=0; ipg<npg; ipg++) {
      if (npg==4) {
	// Four points near  the side centers
	double a=(30.+sqrt(660.))/120.;
	xipg = a; etapg=a; zetapg=a; 
	if (ipg==0) {
	  xipg=1-3*a; 
	} else if (ipg==1) {
	  etapg=1-2*a;
	} else {
	  zetapg=1-2*a;
	}
	wpg[ipg] =  1/12.;

      } else if (npg==1) {
	// One point in the center of the element
	xipg=1./3.;
	etapg=1./3.;
	zetapg=1./3.;
	wpg[ipg] =  1./3.;

      } else {
	PetscPrintf(PETSC_COMM_WORLD,
		    "Not valid value of npg. for tetras %d\n",
		    npg);
	CHKERRA(1); 
      }

      shape[ipg] = RowVector(nel);
      shape[ipg](1)=1-xipg-etapg-zetapg;
      shape[ipg](2)=xipg;
      shape[ipg](3)=etapg;
      shape[ipg](4)=zetapg;
    
      dshapex[ipg]= Matrix(ndimel,nel);
      dshapexi[ipg]= Matrix(ndimel,nel);
    
      dshapexi[ipg]=0;
      dshapexi[ipg](1,2)=1.;
      dshapexi[ipg](2,3)=1.;
      dshapexi[ipg](3,4)=1.;
      dshapexi[ipg](1,1)=-1.;
      dshapexi[ipg](2,1)=-1.;
      dshapexi[ipg](3,1)=-1.;
    }

  } else if (!fnmatch("cartesian*d_face",geom,0)) {

    // This is for integration over the face of an
    // element. (For instance in BCCONV elements.)

    // npg1d:= number of points per dimensional direction
    // npg:= npg1d^ndimel total number of Gauss points
    int ndimel;
    sscanf(geom,"cartesian%dd",&ndimel);

//      if (ndimel==0 && npg==1) {
    int npg1d=npg1d=int(pow(double(npg),1./double(ndimel)));
    // AGREGAR lin1d y brick integrations!!
    //    } else if (!(strcmp(geom,"quad")) || !(strcmp(geom,"lin1d"))
    //  	     || !(strcmp(geom,"brick")) ) {

    if (npg1d!=2) GPERROR;
    // Position of the nodes in the master element
    
    if (ndimel==2) {

      ipg=-1;
      for (int ixipg=0; ixipg<=1; ixipg++) {
	double xipg=(2*ixipg-1)/sqrt(3.);
	double etapg = -1;
	ipg++;
	wpg[ipg] = 1;

	cartesian_2d_shape(shape[ipg],dshapexi[ipg],xipg,etapg);
      }

    } else {
      
      GPERROR;
       
    } 

  } else if (!fnmatch("cartesian*d",geom,0)) {

    // npg1d:= number of points per dimensional direction
    // npg:= npg1d^ndimel total number of Gauss points
    int ndimel;
    sscanf(geom,"cartesian%dd",&ndimel);
    int npg1d=int(pow(double(npg),1./double(ndimel)));
    // AGREGAR lin1d y brick integrations!!
    //    } else if (!(strcmp(geom,"quad")) || !(strcmp(geom,"lin1d"))
    //  	     || !(strcmp(geom,"brick")) ) {

    if (npg1d!=2 && ndimel>0) GPERROR;
    
    if (ndimel==0 && npg==1) {

      shape[0] = RowVector(nel);
      shape[0](1)=1.;
      wpg[0]=1.;
      dshapex[0] = Matrix(ndimel,nel);

    } else if (ndimel==1) {
      // Position of the nodes in the master element
      Matrix xinode(1,2);
      xinode << -1 << 1 ;
      
      ipg=-1;
      for (int ixipg=0; ixipg<=1; ixipg++) {
	double xipg=(2*ixipg-1)/sqrt(3.);
	ipg++;
	wpg[ipg] = 1;
	
	shape[ipg] = RowVector(nel);
	dshapex[ipg]= Matrix(ndimel,nel);
	dshapexi[ipg]= Matrix(ndimel,nel);
	
	for (int iloc=1; iloc<=nel; iloc++) {
	  // 1D Shape Functions
	  // along X
	  double sxi=(1.+xipg*xinode(1,iloc))/2.;
	  double dsxidxi=xinode(1,iloc)/2.;

	  shape[ipg](iloc)=sxi;
	  dshapexi[ipg](1,iloc)=dsxidxi;
	}
      }
    } else if (ndimel==2) {
      Matrix xinode(2,4);
      xinode << -1 << 1 << 1 << -1 << -1 << -1 << 1 << 1;
      
      ipg=-1;
      for (int ixipg=0; ixipg<=1; ixipg++) {
	double xipg=(2*ixipg-1)/sqrt(3.);
	for (int ietapg=0; ietapg<=1; ietapg++) {
	  ipg++;
	  double etapg=(2*ietapg-1)/sqrt(3.);
	  wpg[ipg] = 1;

	  shape[ipg] = RowVector(nel);
	  dshapex[ipg]= Matrix(ndimel,nel);
	  dshapexi[ipg]= Matrix(ndimel,nel);
	
	  for (int iloc=1; iloc<=nel; iloc++) {
	    // 1D Shape Functions
	    // along X
	    double sxi=(1.+xipg*xinode(1,iloc))/2.;
	    double dsxidxi=xinode(1,iloc)/2.;

	    // along Y
	    double seta=(1.+etapg*xinode(2,iloc))/2.;
	    double dsetadeta=xinode(2,iloc)/2.;

	    shape[ipg](iloc)=sxi*seta;
	    dshapexi[ipg](1,iloc)=dsxidxi*seta;
	    dshapexi[ipg](2,iloc)=sxi*dsetadeta;
	  }
	}
      }
    } else if (ndimel==3) {
      Matrix xinode(3,8);
      xinode << -1 << 1 << 1 << -1 << -1 << 1 << 1 << -1 
	     << -1 << -1 << 1 << 1 << -1 << -1 << 1 << 1 
	     << -1 << -1 << -1 << -1 << 1 << 1 << 1 << 1;
	
      
      ipg=-1;
      for (int ixipg=0; ixipg<=1; ixipg++) {
	double xipg=(2*ixipg-1)/sqrt(3.);
	for (int ietapg=0; ietapg<=1; ietapg++) {
	  double etapg=(2*ietapg-1)/sqrt(3.);
	  for (int izetapg=0; izetapg<=1; izetapg++) {
	    double zetapg=(2*izetapg-1)/sqrt(3.);

	    ipg++;
	    wpg[ipg] = 1;

	    shape[ipg] = RowVector(nel);
	    dshapex[ipg]= Matrix(ndimel,nel);
	    dshapexi[ipg]= Matrix(ndimel,nel);
	    
	    for (int iloc=1; iloc<=nel; iloc++) {
	      // 1D Shape Functions
	      // along X
	      double sxi=(1.+xipg*xinode(1,iloc))/2.;
	      double dsxidxi=xinode(1,iloc)/2.;
	      
	      // along Y
	      double seta=(1.+etapg*xinode(2,iloc))/2.;
	      double dsetadeta=xinode(2,iloc)/2.;
	      
	      // along Z
	      double szeta=(1.+zetapg*xinode(3,iloc))/2.;
	      double dszetadzeta=xinode(3,iloc)/2.;
	      
	      shape[ipg](iloc)=sxi*seta*szeta;
	      dshapexi[ipg](1,iloc)=dsxidxi*seta*szeta;
	      dshapexi[ipg](2,iloc)=sxi*dsetadeta*szeta;
	      dshapexi[ipg](3,iloc)=sxi*seta*dszetadzeta;
	    }
	  }
	}
      }
    } else GPERROR;
      
  } else GPERROR;
  
  if (mat_version==GP_FASTMAT) {

    // make fastmat copy version
    FM_shape = new (FastMat *)[npg];
    FM_dshapexi = new (FastMat *)[npg];
    for (int ipg=0; ipg<npg; ipg++) {
      FM_shape[ipg] = new FastMat;
      NM2FM(*FM_shape[ipg],shape[ipg]);
      FM_dshapexi[ipg] = new FastMat;
      NM2FM(*FM_dshapexi[ipg],dshapexi[ipg]);
    }

  } else if (mat_version==GP_FASTMAT2) {

    FM2_shape = new (FastMat2 *)[npg];
    FM2_dshapexi = new (FastMat2 *)[npg];
    for (int ipg=0; ipg<npg; ipg++) {
      FM2_shape[ipg] = new FastMat2;
      FM2_shape[ipg]->set(shape[ipg]);
      int nel = shape[ipg].Ncols();
      int nrows = shape[ipg].Nrows();
      assert(nrows==1);
      FM2_shape[ipg]->reshape(1,nel);
      FM2_dshapexi[ipg] = new FastMat2;
      FM2_dshapexi[ipg]->set(dshapexi[ipg]);
    }

  }
    
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
GPdata::~GPdata() {
  delete[] wpg;
  delete[] shape;
  delete[] dshapexi;
  delete[] dshapex;
  if (mat_version == GP_FASTMAT) {
    for (int ipg=0; ipg<npg; ipg++) {
      delete FM_shape[ipg];
      delete FM_dshapexi[ipg];
      // FM_shape[ipg]->~FastMat();
      // FM_dshapexi[ipg]->~FastMat();
    }
    delete[] FM_shape;
    delete[] FM_dshapexi;
  } else if (mat_version == GP_FASTMAT2) {
    for (int ipg=0; ipg<npg; ipg++) {
      delete FM2_shape[ipg];
      delete FM2_dshapexi[ipg];
      // FM_shape[ipg]->~FastMat();
      // FM_dshapexi[ipg]->~FastMat();
    }
    delete[] FM2_shape;
    delete[] FM2_dshapexi;
  }
}

#undef GPERROR

//__INSERT_LICENSE__
//$Id: mmove.cpp,v 1.2 2002/11/27 19:13:42 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <newmatap.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove.h"

void mesh_move::init() {

  int ierr;
  FastMat2 C;
  double c1 = sqrt(1./3.), c2 = sqrt(1./6.);
  // These are the gradient of shape functions for a master
  // tetra, with unit edge length with nodes at 
  // [+-1/2 0 0], [0,sqrt(3)/2,0], [0,1/sqrt(6),sqrt(2/3)]
  double c[12] = {-1., -c1, -c2, +1., -c1, -c2, 
		  0., 2*c1, -c2, 0., 0., 3*c2};
  
  //o Distortion coefficient
  TGETOPTDEF_ND(thash,double,c_distor,1.);

  //o Functional exponent
  TGETOPTDEF_ND(thash,double,distor_exp,1.);

  //o Distortion coefficient
  TGETOPTDEF_ND(thash,double,c_volume,1.);

  G.resize(2,ndim,ndim); // the metric tensor
  xlocp.resize(1,nel*ndim); // The perturbed coordinates
  xloc0.resize(1,nel*ndim); // The perturbed coordinates
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  assert(nel==ndim+1); // Only for triangles in 2D, tetras in 3D
  J.resize(2,ndim,ndim);
  dNdxi.resize(2,ndim,nel);

  GG.ReSize(ndim);
  D.ReSize(ndim);

  if (ndim==2) {
    dNdxi.setel(-sin(M_PI/3)*cos(M_PI/6),1,1);
    dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,1);
    dNdxi.setel(+sin(M_PI/3)*cos(M_PI/6),1,2);
    dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,2);
    dNdxi.setel(0                       ,1,3);
    dNdxi.setel(+sin(M_PI/3)            ,2,3);
  } else {
#if 0
    dNdxi.set(0.);
    for (int k=1; k<=3; k++) dNdxi.setel(1.,k,k);
    dNdxi.setel(-1.,1,4);
    dNdxi.setel(-1.,2,4);
    dNdxi.setel(-1.,3,4);
#endif
    
    C.resize(2,nel,ndim).set(c).t();
    dNdxi.set(C);
  }
  res_Dir.resize(2,nel,ndim);
  
}

double mesh_move::distor_fun(FastMat2 & xlocp) {
  double df,la1,la2,la3,vol,volref;

  xlocp.reshape(2,nel,ndim);
  J.prod(xlocp,dNdxi,-1,1,2,-1);
  G.prod(J,J,-1,1,-1,2);
  for (int i=1; i<=ndim; i++) {
    for (int j=1; j<=ndim; j++) {
      GG(i,j) = G.get(i,j);
    }
  }
  EigenValues(GG,D);
  la1 = D(1,1);
  la2 = D(2,2);
  if (ndim==3) la3 = D(3,3);

  volref=1.;
  if (ndim==2) {
    vol = la1*la2;
    df = c_distor * pow((la1-la2)*(la1-la2)/(la1*la2),distor_exp)
      + c_volume*(vol-volref)*(vol-volref);
  } else {
    vol = la1*la2*la3;
    double sum = (la1-la2)*(la1-la2) + (la2-la3)*(la2-la3) + (la1-la3)*(la1-la3);
    df = c_distor * pow(sum/pow(vol,2./3.),distor_exp)
      + c_volume*(vol-volref)*(vol-volref);
  }
  xlocp.reshape(1,nel*ndim);

  return df;
}

void mesh_move::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  double distor,distor_p,distor_m,eps=1e-4,
    distor_pp,distor_pm,distor_mp,distor_mm, d2f;
  int nen = nel*ndim;
  xloc0.reshape(2,nel,ndim).set(xloc).reshape(1,nel*ndim);
  distor = distor_fun(xloc0);

  res.reshape(1,nel*ndim);
  for (int k=1; k<=nel*ndim; k++) {
    xlocp.set(xloc0);
    xlocp.setel(xloc0.get(k)+eps,k);
    distor_p = distor_fun(xlocp);

    xlocp.set(xloc0);
    xlocp.setel(xloc0.get(k)-eps,k);
    distor_m = distor_fun(xlocp);

    res.setel(-(distor_p-distor_m)/(2*eps),k);
  }
  res.reshape(2,nel,ndim);

  mat.reshape(2,nen,nen);
  for (int k=1; k<=nel*ndim; k++) {
    for (int l=k; l<=nel*ndim; l++) {
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)+eps,k);
      xlocp.setel(xlocp.get(l)+eps,l);
      distor_pp = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)+eps,k);
      xlocp.setel(xlocp.get(l)-eps,l);
      distor_pm = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)-eps,k);
      xlocp.setel(xlocp.get(l)+eps,l);
      distor_mp = distor_fun(xlocp);
      
      xlocp.set(xloc0);
      xlocp.setel(xloc0.get(k)-eps,k);
      xlocp.setel(xlocp.get(l)-eps,l);
      distor_mm = distor_fun(xlocp);

      d2f = (distor_pp - distor_pm - distor_mp + distor_mm)/(4.*eps*eps);
      mat.setel(d2f,k,l);
    }
  }

  for (int k=2; k<=nel*ndim; k++) {
    for (int l=1; l<=k-1; l++) {
      d2f = mat.get(l,k);
      mat.setel(d2f,k,l);
    }
  }

  mat.reshape(4,nel,ndim,nel,ndim);
  res_Dir.prod(mat,state_new,1,2,-1,-2,-1,-2);
  res.axpy(res_Dir,-1.);
    
}

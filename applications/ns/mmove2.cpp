//__INSERT_LICENSE__
//$Id: mmove2.cpp,v 1.4 2002/12/03 23:35:02 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include <src/texthf.h>
#ifdef USE_NEWMAT
#include <newmatap.h>
#endif

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove.h"

extern GlobParam *GLOB_PARAM;

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::init() {
  int ierr;
  FastMat2 C;
  double c1 = sqrt(1./3.), c2 = sqrt(1./6.);
  // These are the gradient of shape functions for a master
  // tetra, with unit edge length with nodes at 
  // [+-1/2 0 0], [0,sqrt(3)/2,0], [0,1/sqrt(6),sqrt(2/3)]
  double c[12] = {-1., -c1, -c2, +1., -c1, -c2, 
		  0., 2*c1, -c2, 0., 0., 3*c2};
  
  G.resize(2,ndim,ndim); // the metric tensor
  xlocp.resize(1,nel*ndim); // The perturbed coordinates
  xloc0.resize(1,nel*ndim); // The perturbed coordinates
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  // assert(nel==ndim+1); // Only for triangles in 2D, tetras in 3D
  J.resize(2,ndim,ndim);
  dNdxi.resize(2,ndim,nel);

  if (ndim==2) {
    if (nel==3) {
      dNdxi.setel(-sin(M_PI/3)*cos(M_PI/6),1,1);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,1);
      dNdxi.setel(+sin(M_PI/3)*cos(M_PI/6),1,2);
      dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,2);
      dNdxi.setel(0                       ,1,3);
      dNdxi.setel(+sin(M_PI/3)            ,2,3);
    } else if (nel==4) {
      double cq[] = {-1,-1,1,-1,1,1,-1,1};
      C.resize(2,nel,ndim).set(cq).t();
      dNdxi.set(C);
    } else PETSCFEM_ERROR("Only tringles ad quads in 2D: nel %d\n",nel);
  } else {
    C.resize(2,nel,ndim).set(c).t();
    dNdxi.set(C);
  }
  glambda.resize(3,ndim,nel,ndim);
  dFdl.resize(1,ndim);
  d2Fdl2.resize(2,ndim,ndim);
  eps = 1e-5;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::la_grad(const FastMat2 &x,FastMat2 &lambda,
				 FastMat2 &glambda) {
  J.prod(x,dNdxi,-1,1,2,-1);
  G.prod(J,J,-1,1,-1,2);
  lambda.seig(G,V);
  tmp1.prod(J,V,1,-1,-1,2);
  tmp2.prod(dNdxi,V,-1,1,-1,2);
  for (int q=1; q<=ndim; q++) {
    glambda.ir(1,q);
    tmp1.ir(2,q);
    tmp2.ir(2,q);
    glambda.prod(tmp1,tmp2,2,1);
  }
  glambda.rs().scale(2.);
  tmp1.rs();
  tmp2.rs();
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
double mesh_move_eig_anal::dfun(const FastMat2 &D) {
  double F=0;
  double vol=1.;
  for (int k=1; k<=ndim; k++) {
    vol *= D.get(k);
    for (int l=1; l<=ndim; l++) 
      F += square(D.get(k)-D.get(l));
  }
  F /= pow(vol,2./double(ndim));
  return F;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::df_grad(const FastMat2 &x,FastMat2 &dFdx) {
  la_grad(x,lambda,glambda);
  double F, F0 = dfun(lambda);
  for (int k=1; k<=ndim; k++) {
    lambdap.set(lambda).addel(eps,k);
    F  = dfun(lambdap);
    dFdl.setel((F-F0)/eps,k);
  }
  dFdx.prod(dFdl,glambda,-1,-1,1,2);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::element_connector(const FastMat2 &xloc,
					   const FastMat2 &state_old,
					   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  const FastMat2 *state = (glob_param->inwt ? &state_new : &state_old);
  x0.set(xloc).add(*state);
  df_grad(x0,res);
  x0.reshape(1,nel*ndim);
  mat.reshape(3,nel,ndim,nel*ndim);
  for (int k=1; k<=nel*ndim; k++) {
    xp.set(x0).addel(eps,k).reshape(2,nel,ndim);
    df_grad(xp,resp);
    xp.reshape(1,nel*ndim);
    mat.ir(3,k).set(resp).rest(res).scale(1/eps);
  }
  x0.reshape(2,nel,ndim);
  mat.rs().reshape(4,nel,ndim,nel,ndim);

  if (!glob_param->inwt) {
    x0.set(xloc).add(state_new);
    df_grad(x0,res);
  } else {
    double q=1;
  }
  res.scale(-1.);
}

//__INSERT_LICENSE__
//$Id: mmove2.cpp,v 1.2 2002/12/02 01:00:32 mstorti Exp $

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
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>---: 
void mesh_move_eig_anal::df_grad(const FastMat2 &x,FastMat2 &lambda,
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
void mesh_move_eig_anal::element_connector(const FastMat2 &xloc,
					   const FastMat2 &state_old,
					   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  df_grad(xloc,lambda,glambda);
  x0.set(xloc);
  x0.reshape(1,nel*ndim);
  double epsilon=1e-4;
  glambda_diff.set(glambda).set(0.).reshape(2,ndim,nel*ndim);
  glp.set(glambda);
  for (int k=1; k<=nel*ndim; k++) {
    xp.set(x0);
    xp.ir(1,k).add(epsilon).rs();
    xp.reshape(2,nel,ndim);
    df_grad(xp,lambdap,glp);
    glambda_diff.ir(2,k).set(lambdap).rest(lambda).scale(1/epsilon);
    xp.reshape(1,nel*ndim);
  }
  x0.reshape(2,nel,ndim);
  glambda_diff.rs().reshape(3,ndim,nel,ndim);
}

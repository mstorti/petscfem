//__INSERT_LICENSE__
//$Id: mmove.cpp,v 1.1.2.1 2001/12/20 02:32:24 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove.h"

void mesh_move::init() {

  int ierr;
  G.resize(1,ndim); // the metric tensor
  xlocp.resize(2,nel,ndim); // The perturbed coordinates
  assert(ndof==ndim);
  assert(ndim==2 || ndim==3);
  assert(nel==ndim+1); // Only for triangles in 2D, tetras in 3D
  J.resize(2,ndim,ndim);
  dNdxi.resize(2,ndim,nel);
  assert(ndim==2);

  dNdxi.setel(-sin(M_PI/3)*cos(M_PI/6),1,1);
  dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,1);
  dNdxi.setel(+sin(M_PI/3)*cos(M_PI/6),1,2);
  dNdxi.setel(-sin(M_PI/3)*sin(M_PI/6),2,2);
  dNdxi.setel(0                       ,1,3);
  dNdxi.setel(+sin(M_PI/3)            ,2,3);
  
}

double mesh_move::distor_fun(FastMat2 & xlocp) {
  double df;
  J.prod(xlocp,dNdxi,-1,1,2,-1);
  G.prod(J,J,-1,1,-1,2);
  return df;
}

void mesh_move::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){

  double distor;
  xlocp.set(xloc);
  distor = distor_fun(xlocp);
    
}

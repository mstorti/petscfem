//__INSERT_LICENSE__
//$Id: mmove2.cpp,v 1.15 2005/06/18 12:02:30 mstorti Exp $

#include <src/fem.h>
#include <src/utils.h>
#include <src/readmesh.h>
#include <src/getprop.h>
#include <src/fastmat2.h>
#include "./fm2funm.h"

#include "nsi_tet.h"
#include "adaptor.h"
#include "mmove2.h"

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mesh_move2::init() {
  assert(ndof==ndim);
  assert(nel==ndim+1);
  x.resize(2,nel,ndim);
  nedge = nel*(nel-1)/2;
  xedge.resize(1,ndim);
  A.resize(2,nel,ndim+1);
  A.ir(2,ndim+1).set(1.0).rs();
  vol_coef = (ndim==3 ? 1./6. : 1./2.);
  qmax = sqrt(3.0)/4.0;
  assert(ndim==2); // not computed yet best quality for ndim==3
  nen = nel*ndim;
  g.resize(1,nen);
  gp.resize(1,nen);
  gm.resize(1,nen);
  xp.resize(1,nen);
  epsi = 1e-6;
  chard = 1;
  qmin = 0.3;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double mesh_move2::hardfun(double q) {
  double coef = (chard-1.0)/(chard+1.0);
  double r = q/qmin;
  double w;
  if (fabs(r)<1e-7) w = 1.0/(1.0-r*r/3.0);
  else w = r/tanh(r);
  return -q + coef*qmin*w;
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
double mesh_move2::disfun(FastMat2 &x) {
  int edge = 0;
  double sum_ledge2=0.;
  x.reshape(2,nel,ndim);
  for (int k=1; k<nel; k++) {
    for (int l=k+1; l<=nel; l++) {
      edge++;
      x.ir(1,k);
      xedge.set(x);
      x.ir(1,l);
      xedge.rest(x);
      sum_ledge2 += xedge.sum_square_all();
    }
  }
  double ledge = sqrt(sum_ledge2/nedge);
  x.rs();
  A.is(2,1,ndim).set(x).rs();
  double vol = A.det()*vol_coef;
  double q = vol/pow(ledge,ndim)/qmax;
  x.reshape(1,nen);
  return hardfun(q);
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mesh_move2::gdisfun(FastMat2 &x, FastMat2 &g) {
  x.reshape(1,nel*ndim);
  
  for (int j=1; j<=nen; j++) {
    xp.set(x).addel(epsi,j);
    double fp = disfun(xp);
    xp.set(x).addel(-2*epsi,j);
    double fm = disfun(xp);
    g.setel((fp-fm)/(2*epsi),j);
  }
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mesh_move2::element_connector(const FastMat2 &xloc,
				   const FastMat2 &state_old,
				   const FastMat2 &state_new,
				   FastMat2 &res,FastMat2 &mat){
  x.set(xloc).add(state_new);
  gdisfun(x,g);
  g.reshape(2,nel,ndim);
  res.set(g).scale(-1.0);
  
  g.reshape(1,nen);
  g.print("");
  mat.reshape(2,nen,nen);
  for (int j=1; j<=nen; j++) {
    xp.set(x).addel(epsi,j);
    gdisfun(xp,gp);

    xp.set(x).addel(-2*epsi,j);
    gdisfun(xp,gm);

    mat.ir(2,j).set(gp).rest(gm).scale(1./(2*epsi));
  }
  mat.rs().reshape(4,nel,ndof,nel,ndof);
#if 1
  mat.print(6,"");
#endif
}

//__INSERT_LICENSE__
//$Id: mmove2.cpp,v 1.16 2005/06/21 01:36:16 mstorti Exp $

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
  int ierr;
  //o Hardening coefficient for inverted element. 
  TGETOPTDEF_ND(thash,double,chard,10.);
  assert(ndof==ndim);
  assert(nel==ndim+1);
  x.resize(2,nel,ndim);
  nedge = nel*(nel-1)/2;
  xedge.resize(1,ndim);
  dx.resize(1,ndim);
  s.resize(1,ndim);
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
#if 0
  chard = 1;
  qmin = 0.3;
#else
  qmin = 0.3;
#endif
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
#if 0
  int edge = 0;
  double sum_ledge2=0.;
  x.reshape(2,nel,ndim);
  for (int k=1; k<nel; k++) {
    for (int l=k+1; l<=nel; l++) {
      edge++;
      x.ir(1,k);
      xedge.set(x);
      x.ir(1,l);
      xedge.minus(x);
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
#elif 0
  assert(ndim==2);
  double sum_l2=0.;
  double l2max = 0., f=0.;
  x.reshape(2,nel,ndim);
  for (int k=0; k<nel; k++) {
    int n1 = k+1;
    int n2 = (k+1)%nel+1;
    int n3 = (k+2)%nel+1;
    x.ir(1,n2);
    xedge.set(x);
    x.ir(1,n3);
    dx.set(x);
    x.ir(1,n1);
    xedge.minus(x);
    dx.minus(x);
    double l2 = xedge.sum_square_all();
    sum_l2 += l2;
    s.setel(-xedge.get(2),1);
    s.setel(+xedge.get(1),2);
    tmp.prod(dx,s,-1,-1);
    double l = sqrt(l2);
    double ff = double(tmp)/l-qmin*l;
    if (l2>l2max) {
      l2max = l2;
      f = chard*ff*ff*(ff<0);
    }
  }
  x.reshape(1,nen);
  f += sum_l2;
  return f;
#else
  assert(ndim==2);
  double l2max = 0., f=0.;
  x.reshape(2,nel,ndim);
  for (int k=0; k<nel; k++) {
    int n1 = k+1;
    int n2 = (k+1)%nel+1;
    int n3 = (k+2)%nel+1;
    x.ir(1,n2);
    xedge.set(x);
    x.ir(1,n3);
    dx.set(x);
    x.ir(1,n1);
    xedge.minus(x);
    dx.minus(x);
    double l2 = xedge.sum_square_all();
    s.setel(-xedge.get(2),1);
    s.setel(+xedge.get(1),2);
    tmp.prod(dx,s,-1,-1);
    double l = sqrt(l2);
    double ff = double(tmp)/l-qmin*l;
    f += l2 + chard*ff*ff*(ff<0);
  }
  x.reshape(1,nen);
  return f;
#endif
}

//---:---<*>---:---<*>---:---<*>---:---<*>---:---<*>
void mesh_move2::gdisfun(FastMat2 &x, FastMat2 &g) {
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

#if 1
  double xx[] = {-1,0,1,0,0,0};
  x.reshape(1,nel*ndim).set(xx);
  double L = 2.;
  int N=41;
  for (int j=0; j<=N; j++) {
    double X = -L+j*2*L/N;
    for (int k=0; k<=N; k++) {
      double Y = -L+k*2*L/N;
      x.setel(X,5).setel(Y,6);
      printf("%g %g %g\n",X,Y,disfun(x));
    }
  }
  PetscFinalize();
  exit(0);
#endif  
  
  x.reshape(1,nel*ndim);
  gdisfun(x,g);
  g.reshape(2,nel,ndim);
  res.set(g).scale(-1.0);

  g.reshape(1,nen);
  // g.print("");
  mat.reshape(2,nen,nen);
  for (int j=1; j<=nen; j++) {
    xp.set(x).addel(epsi,j);
    gdisfun(xp,gp);

    xp.set(x).addel(-2*epsi,j);
    gdisfun(xp,gm);

    mat.ir(2,j).set(gp).minus(gm).scale(1./(2*epsi));
  }
  mat.rs().reshape(4,nel,ndof,nel,ndof);
#if 0
  mat.print(6,"");
#endif
}

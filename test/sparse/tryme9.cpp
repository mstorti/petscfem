/*__INSERT_LICENSE__*/
// $Id: tryme9.cpp,v 1.5 2003/02/28 23:51:05 mstorti Exp $

#include <src/cloud2.h>

int main(int argc, char **argv) {
  Cloud2 cloud;
  int nderiv = 7, nx = 9, ndim=2;
  int derivs[] = {0,0,0,1,1,0,0,2,2,0,1,1,2,2};
  int npol[] = {2,2};
  FastMat2 x(3,3,3,2), x0(1,2), w(2,nx,nderiv);
  for (int j=1; j<=3; j++) {
    for (int k=1; k<=3; k++) {
      x.setel(double(j)-2,j,k,1);
      x.setel(double(k)-2,j,k,2);
    }
  }
  x.reshape(2,nx,ndim);
  cloud.init(ndim,nx,nderiv,derivs,npol);
  cloud.coef(x,w);
  w.print();

  x0.set(1.);
  x.add(1.);
  cloud.coef(x,w,x0);
  w.print();

}

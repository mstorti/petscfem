/*__INSERT_LICENSE__*/
// $Id: tryme9.cpp,v 1.1 2003/02/27 03:32:41 mstorti Exp $

#include <src/cloud2.h>

int main(int argc, char **argv) {
  Cloud2 cloud;
  int derivs[] = {1,0,0,1,0,0};
  int npol[] = {2,2};
  FastMat2 x(3,3,3,2), x0(1,2), w(2,9,3);
  for (int j=1; j<=3; j++) {
    for (int k=1; k<=3; k++) {
      x.setel(double(j)-2,j,k,1);
      x.setel(double(k)-2,j,k,2);
    }
  }
  x.reshape(2,9,2);
  cloud.init(2,9,3,derivs,npol);
  cloud.coef(x,w,x0);
  w.print("w: ");

}

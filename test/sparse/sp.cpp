/*__INSERT_LICENSE__*/
// $Id: sp.cpp,v 1.3 2001/09/21 15:37:57 mstorti Exp $

#include <cmath>
#include <vector>

#include <sparse.h>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

using namespace Sparse;

#define M 1000
int main() {
  int j,k,m;

  Vec v(5),w(4),u;
  Mat a;
  Indx I(14,17),J(8,14,2);

  v.set(2,1.).set(3,2.).print("set values ...  ");
  v.set(10,10.).set(11,11.).print("set at 10 11 ...  ");
  v.grow(0).set(5,5.).print("set at 5 ...  ");
  // v.set(20,20.).print(); should give an error
  v.resize(10).print("resized to 10... ");
  v.clear().print("cleared ...  ");
  w.set(0,14.).set(1,15.).set(2,16.).set(3,17.).print("w: ");
  v.resize(20).set(I,w).print("setting v[I] = w  ");

  u.resize(20).set(J,v,I).print("u[I] = v[J]");
  u.scale(2.).print("u*2 = ");

  u.print("u:");
  v.print("v:");
  u.axpy(3.,v).print("u = u+3*v ... ");

  m=5;
  a.resize(m,m);
  for (j=0; j<m; j++) 
    for (k=0; k<m; k++) 
      if ((j+k) % 2) a.set(j,k,double(j*10+k));
  a.print("checkerboard mat");

}

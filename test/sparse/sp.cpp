/*__INSERT_LICENSE__*/
// $Id: sp.cpp,v 1.6 2001/09/21 22:47:40 mstorti Exp $

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

void my_set(Vec &v,int k=1) {
  v.clear();
  int l = v.length();
  for (int j=0; j<l; j+= k) 
    v.set(j,double(j));
}

#define M 1000
int main() {
  int j,k,m;

  Vec v(5),w(4),u;
  Mat a,b;
  Indx I(14,17),J(8,14,2),K(0,4,2);

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

  printf("u empty? %d\n",u.empty());
  u.clear();
  printf("u empty after clear() ? %d\n",u.empty());

  m=5;
  a.resize(m,m);
  for (j=0; j<m; j++) 
    for (k=0; k<m; k++) 
      if ((j+k) % 2) a.set(j,k,double(j*10+k));
  a.print_f("checkerboard mat");
  
  u.resize(5);
  my_set(u);
  u.print("after my_set; ");
  for (j=0; j<2; j++) 
    a.setr(j,u);
  a.print_f("rows 0 1  set to u: ");

  for (j=3; j<5; j++) 
    a.setc(j,u);
  a.print_f("cols 3 4 set to u: ");

  u.set(a,3).print("u = 3rd row of a\n");
  u.setc(a,3).print("u = 3rd col of a\n");

  b.setr(a,K).print_f("b set to rows 0 2 4 of a");
  b.setc(a,K).print_f("b set to cols 0 2 4 of a");

}

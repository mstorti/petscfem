//__INSERT_LICENSE__
// $Id: project.cpp,v 1.2 2005/02/24 01:20:03 mstorti Exp $

#include <cstdio>
#include <src/fastmat2.h>
#include <src/dvector.h>
#include <src/dvector2.h>

int main() {
  int ndim = 3;
  int nel = 3;

  dvector<double> x1,x2;
  dvector<int> ico1;

#define DATA_DIR "/home/mstorti/PETSC/tatuus/data/"

  x1.cat(DATA_DIR "static_p_blade.nod").defrag();
  assert(x1.size() % ndim ==0);
  int nnod1 = x1.size()/ndim;
  x1.reshape(2,nnod1,ndim);

  ico1.cat(DATA_DIR "blade.con").defrag();
  assert(ico1.size() % nel ==0);
  int nelem1 = ico1.size()/nel;
  ico1.reshape(2,nelem1,nel);

  printf("%d nodes, %d elems read\n",nnod1,nelem1);

  dvector<double> Sinv_v;
  Sinv_v.a_resize(3,nelem1,ndim,ndim).defrag();

  FastMat2 S(2,ndim,ndim),Sinv(2,ndim,ndim),xn1(1,ndim),
    xn2(1,ndim),xn3(1,ndim),x12(1,ndim),x13(1,ndim);
  for (int k=0; k<nelem1; k++) {
    // Nodes of element (base 0)
    int 
      n1 = ico1.e(k,0),
      n2 = ico1.e(k,1),
      n3 = ico1.e(k,2);
    // Position of first node
    xn1.set(&x1.e(n1-1,0)) ;
    xn2.set(&x1.e(n2-1,0)) ;
    xn3.set(&x1.e(n3-1,0)) ;
    x12.set(xn2).rest(xn1);
    S.ir(1,1).set(x12);
    x13.set(xn3).rest(xn1);
    S.ir(1,2).set(x13);
    S.ir(1,3).cross(x12,x13).rs();
    Sinv.inv(S);
    Sinv.export_vals(&Sinv_v.e(k,0,0));
  }

  for (int 

}

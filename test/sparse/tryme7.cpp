/*__INSERT_LICENSE__*/
// $Id: tryme7.cpp,v 1.2 2003/02/17 01:28:02 mstorti Exp $

#include <src/dvector.h>
#include <src/dvector2.h>
#include <cstdio>

void dvprint(dvector<int> &v,int M,int N) {
  for (int j=0; j<M; j++) {
    printf("v(%d,*) = ",j);
    for (int k=0; k<N; k++) printf(" %d",v.e(j,k));
    printf("\n");
  }
}

int main(int argc, char **argv) {
  dvector<int> v;
  int M=10,N=10;
  v.a_resize(2,M,N);
  for (int j=0; j<v.size(); j++) v.ref(j)=j;
  dvprint(v,M,N);

  printf("\n\n");
  M=20; N=5;
  v.reshape(2,M,N);
  dvprint(v,M,N);
  
}

/*__INSERT_LICENSE__*/
// $Id: tryme8.cpp,v 1.1 2003/02/25 20:34:25 mstorti Exp $

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
  dvector<int> v(10);
  int M=10,N=10;
  v.a_resize(2,M,N);
  for (int j=0; j<v.size(); j++) v.ref(j)=j;
  dvprint(v,M,N);
  printf("number f chunks: %d\n",v.chunks_n());

  // v.set_chunk_size(v.size());//defrag
  v.defrag();
  printf("After defrag\n");
  dvprint(v,M,N);
  printf("number f chunks: %d\n",v.chunks_n());
  
}

/*__INSERT_LICENSE__*/
// $Id: tryme8.cpp,v 1.2 2003/02/27 03:32:41 mstorti Exp $

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

  // export values to a C array, scale by 10 and
  // then get back to the dvector
  int *array = new int[M*N];
  v.export_vals(array);
  for (int j=0; j<M*N; j++) array[j] *= 10;
  v.set(array);
  printf("After mult by 10\n");
  dvprint(v,M,N);

  // v.set_chunk_size(v.size());//defrag
  v.defrag();
  printf("After defrag\n");
  dvprint(v,M,N);
  printf("number f chunks: %d\n",v.chunks_n());

  delete[] array;
}

/*__INSERT_LICENSE__*/
// $Id: tryme7.cpp,v 1.1 2003/02/16 22:03:18 mstorti Exp $

#include <src/dvector.h>
#include <cstdio>

int main(int argc, char **argv) {
  dvector<int> v;
  v.a_resize(2,2,2);
  for (int j=0; j<v.size(); j++) v.ref(j)=j;
  for (int j=0; j<2; j++) {
    printf("v(%d,*) = ",j);
    for (int k=0; k<2; k++) printf(" %d",v.e(j,k));
    printf("\n");
  }
}

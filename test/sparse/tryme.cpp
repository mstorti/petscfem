/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.1 2001/09/21 15:09:32 mstorti Exp $

// Example of use of `lower_bound'. 

#include <cmath>
#include <cstdio>

#include <map>
//#include <algorithm>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

#define M 100
#define N 1000

#if 0
class Trunc {
public:
  int operator()(pair<int,double> p) {return p.first < N/2;}
};
#endif

int main() {
  int j,k;
  double a;
  map<int,double> m;
  map<int,double>::iterator i;

  for (j=0; j<M; j++) {
    a = drand();
    k = irand(1,N);
    m[k] = a;
  }

  printf("max val %d, number of elements %d\nAll elements: \n",N,M);
  for (i=m.begin(); i!=m.end(); i++) 
    printf("%d -> %f\n",i->first,i->second);

  i = m.lower_bound(N/2);
  m.erase(i,m.end());

  printf("=============\n\nTruncate at %d \n",N/2);
  for (i=m.begin(); i!=m.end(); i++) 
    printf("%d -> %f\n",i->first,i->second);
}

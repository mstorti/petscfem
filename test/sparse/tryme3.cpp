/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.2 2002/07/18 18:33:12 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <src/utils.h>
#include <src/iisdgraph.h>

TextHashTable GLOBAL_OPTIONS;

#if 0
double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}
#endif

int main() {
  StoreGraph g;
  int M=1000000;
  for (int j=0; j<M; j++) { g.add(irand(1,M*M),irand(1,M*M)); }
  g.clear();
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  delete[] a;
}

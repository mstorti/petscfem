/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.1 2002/07/18 18:10:15 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <src/graph>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

int main() {
  StoreGraph g;
  int M=1000000;
  for (int j=0; j<M; j++) { g.add(irand(M*M),irand(M*M)); }
  s.clear();
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  delete[] a;
}

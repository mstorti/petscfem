/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.5 2002/07/18 22:05:41 mstorti Exp $

#include <cmath>
#include <map>
#include <set>
#include <vector>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

//#define ALLOC malloc_alloc
#define ALLOC alloc
typedef set< int, less<int>, ALLOC > SET;

class MAP : public map< int, SET, less<int>, ALLOC > {
public:
  void add(int i,int j) { (*this)[i].insert(j); }
};

int main(int argc, char **argv) {
  int M=1000000;
  MAP g;
  for (int j=0; j<M; j++) { g.add(irand(1,M),irand(1,M)); }
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  g.clear();
  delete[] a;
}

/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.4 2002/07/21 14:50:08 mstorti Exp $

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

int main() {
  int jjjj=1;
#if 1
  int M=300000;
  typedef set<double,less<double>,malloc_alloc> S;
  // typedef set<double> S;
  S s;
  double *a = new double[3*M];
  int ia=0;
  for (int j=0; j<M; j++) { 
    s.insert(drand());
    for (int jj=0; jj<3; jj++) a[ia++] = double(ia);
  }
  s.clear();
  delete[] a;
#elif 0
  int M=1000000;
  vector<double> s;
  s.resize(3*M);
  for (int j=0; j<3*M; j++) s[j] = double(j);
  s.clear();
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  delete[] a;
#else
  int M=3000000;
  map<int,double> s;
  for (int j=0; j<M; j++) {
    int j = 100*irand(1,M);
    s.insert(pair<int,double>(j,drand()));
  }
  s.clear();
  double *a = new double[2*M];
  for (int j=0; j<2*M; j++) a[j] = double(j);
  delete[] a;
  a=NULL;
#endif
}

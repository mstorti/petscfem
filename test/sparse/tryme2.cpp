/*__INSERT_LICENSE__*/
// $Id: tryme2.cpp,v 1.1 2002/07/17 20:11:49 mstorti Exp $

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
#if 0
  int M=1000000;
  set<double> s;
  for (int j=0; j<M; j++) { s.insert(drand()); }
  s.clear();
  double *a = new double[3*M];
  for (int j=0; j<3*M; j++) a[j] = double(j);
  delete[] a;
#elif 1
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

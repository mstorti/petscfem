//__INSERT_LICENSE__
// $Id: tryme2.cpp,v 1.1 2004/12/28 14:06:51 mstorti Exp $

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

int its;

int ilog2(int x) {
  return int(ceil(log(double(x))/log(2.0)));
}

int lower_bound(vector<int> &v, int x) {
  its = 0;
  int 
    j, vj,
    j1 = 0,
    v1 = v[j1],
    j2 = v.size()-1,
    v2 = v[j2];
  if (x > v2) return v.size();
  if (x == v2) {
    while (j2>=0 && v[j2]==x) { j2--; }
    return j2;
  }
  while (j2-j1 > 1) {
    its++;
    if (v2 == v1) return j1;
    j = j1 + int((double(x-v1)*double(j2-j1))/double(v2-v1));
    if (j==j1) j++;
    if (j==j2) j--;
    vj = v[j];
    if (x < vj ) { j2=j; v2=vj; }
    else if (x > vj ) { j1=j; v1=vj; }
    else {
      while (j>0 && v[j-1]==x) { j--; }
      return j;
    }
  }
}

int lower_bound2(vector<int> &v, int x) {
  its = 0;
  int 
    j, 
    j1 = 0,
    j2 = v.size();
  while (j2-j1 > 1) {
    its++;
    j = j1 + (j2-j1)/2;
    vj = v[j];
    if (x < vj ) { j2 = j; }
    else if (x >= vj ) { j1 = j; }
    else {
      while (j>0 && v[j-1]==x) { j--; }
      return j;
    }
  }
}

int main() {
  int N=1000000;
  vector<int> v(N);
  for (int j=0; j<N; j++)
    v[j] = rand();
  sort(v.begin(),v.end());

  int its=
  for (int k=0; k<10; k++) {
    int j = rand() % N;
    int jj = lower_bound(v,v[j]);
    assert(v[jj]==v[j]);
    printf("%d its, bin %d\n",its,ilog2(N));
  }
}

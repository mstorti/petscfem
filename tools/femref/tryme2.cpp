//__INSERT_LICENSE__
// $Id: tryme2.cpp,v 1.3 2004/12/28 23:27:02 mstorti Exp $

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

#if 0
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
#endif

static int 
lower_boundb_aux(vector<int> &v, 
	     int x, int j) {
  while (j>=0 && v[j]==x) { j--; }
  return j;
}

int lower_boundb(vector<int> &v, 
		 int x, int j1, int j2) {
  while (j2-j1 > 1) {
    its++;
    int j = j1 + (j2-j1)/2;
    int vj = v[j];
    if (x <= vj ) j2 = j;
    else  j1 = j; 
  }
  return j2;
}

int lower_boundbl(vector<int> &v, int x) {
  its = 0;
  int 
    j, vj,
    j1 = 0,
    v1 = v[j1],
    j2 = v.size()-1,
    v2 = v[j2];
  if (x > v2) return v.size();
  if (x == v2 || j2-j1 <= 1) 
    return lower_boundb_aux(v,x,j2);
  j = j1 + int((double(x-v1)*double(j2-j1))/double(v2-v1));
  int dj = int(sqrt(double(j)));
  int jj1 = j-dj;
  if (jj1>=0) {
    int vjj1 = v[jj1];
    if (vjj1 <= x) {
      v1=vjj1;
      j1=jj1;
    }
  }
  int jj2 = j + dj;
  if (jj2<v.size()) {
    int vjj2 = v[jj2];
    if (vjj2 > x) {
      v2 = vjj2;
      j2=jj2;
    }
  }
  return lower_boundb(v,x,j1,j2);
}

int lower_boundb(vector<int> &v, int x) {
  its = 0;
  if (v[0]>x) return 0;
  return lower_boundb(v,x,0,v.size());
}

int main() {
  int N=10000;
  vector<int> v(N);
  for (int j=0; j<N; j++)
    v[j] = rand();
  sort(v.begin(),v.end());

  for (int k=0; k<10; k++) {
    int j = rand() % N;
    // int jj = lower_boundb(v,v[j]);
    int jj = lower_boundbl(v,v[j]);
    assert(v[jj]==v[j]);
    printf("%d its, bin %d\n",its,ilog2(N));
  }
}

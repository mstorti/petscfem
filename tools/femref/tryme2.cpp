//__INSERT_LICENSE__
// $Id: tryme2.cpp,v 1.7 2004/12/29 00:52:30 mstorti Exp $

#include <ctime>
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

static int 
lower_boundb_aux(vector<int> &v, 
	     int x, int j) {
  while (j>=0 && v[j]==x) { j--; }
  return j;
}

int lower_boundb(vector<int> &v, 
		 int x, int j1, int j2) {
  // printf("j1 %d, j2 %d\n",j1,j2);
  while (j2-j1 > 1) {
    its++;
    int j = j1 + (j2-j1)/2;
    int vj = v[j];
    if (x <= vj ) j2 = j;
    else  j1 = j; 
  }
  return j2;
}

int DJGLOB;

int lower_boundbl(vector<int> &v, int x) {
  its = 0;
  int 
    j, vj,
    j1 = 0,
    v1 = v[j1],
    j2 = v.size()-1,
    v2 = v[j2];
  if (x > v2) return v.size();
  if (x < v1) return 0;
  if (v1 == v2) return 0;
  j = j1 + int((double(x-v1)*double(j2-j1))/double(v2-v1));
  // int dj = int(pow(double(j),0.5));
  int dj = int(sqrt(double(j)));
  // int dj = DJGLOB;
  int jj1 = j-dj;
  if (jj1>=0 && v[jj1] < x) j1=jj1;
  int jj2 = j + dj;
  if (jj2<v.size() && v[jj2]>= x) j2=jj2;
  return lower_boundb(v,x,j1,j2);
}

int lower_boundb(vector<int> &v, int x) {
  its = 0;
  if (v[0]>=x) return 0;
  // printf("size %d\n",v.size());
  return lower_boundb(v,x,0,v.size());
}

int main() {
  time_t start, end;
  int N=10000000, ntries=10000000;
  DJGLOB = int(sqrt(double(N)));
  vector<int> v(N);
  for (int j=0; j<N; j++) 
    v[j] = rand();
  sort(v.begin(),v.end());

#if 0
  for (int j=0; j<1000; j++) 
    printf("v[%d]=%d\n",j,v[j]);
#endif

  char *method = NULL;
#define USE_LINEAR

#ifndef USE_LINEAR
    method = "binary search";
#define FUN lower_boundb
#else
    method = "linear pred. + binary search";
#define FUN lower_boundbl
#endif
  
  start = time(NULL);
  int iters = 0;
  for (int k=0; k<ntries; k++) {
    int j = rand() % N;

    int jj = FUN(v,v[j]);
    iters += its;
    if (v[jj]!=v[j]) {
      printf("v[j=%d]=%d, v[jj=%d]=%d, k %d\n",
	     j,v[j],jj,v[jj],k);
    }
    // printf("%d its, bin %d\n",its,ilog2(N));
  }
  double elaps = difftime(time(NULL),start);
  printf("Stats for %s method:\n" 
	 "%d tries, %d elements\n"
	 "%d its, rate %fits/try\n"
	 "%d secs, rate %g secs/T-elem-try\n",
	 method, ntries, N, iters, 
	 double(iters)/double(ntries),
	 int(elaps), 
	 elaps/(double(N)*double(ntries))*1e12);

}

//__INSERT_LICENSE__
// $Id: trymem.cpp,v 1.2 2004/12/26 00:41:00 mstorti Exp $

#include <cstdio>
#include <cassert>
#include <sys/resource.h>
#include <sys/vtimes.h>
#include <sys/unistd.h>
#include <cstdlib>
#include <vector>

using namespace std;

struct A {
  int k;
  double av[2];
};


void report(int N,int Asize) {
  char file[1000], *line;
  int  mem;
  size_t n=0;
  sprintf(file,"/proc/%d/status",getpid());
  FILE *fid = fopen(file,"r");
  while(1) {
    assert(getline(&line,&n,fid)!=-1);
    if (sscanf(line,"VmRSS: %d kB",&mem)) break;
  }
  fclose(fid);
  double Aused = double(mem)*1024.0/N;
  printf("%dK used, %f bytes/struct approx, excedent %f\n",mem,
	 Aused,Aused-Asize);
}

int main() {
  int N = 10000000;
  int Asize = sizeof(A);
  vector<A> va;
  A *ap;
#if 1
  va.resize(N);
  printf("------------------\nvector<A>:\n");
#endif
#if 0
  printf("------------------\nnew A (A cells):\n");
  for (int j=0; j<N; j++) new A;
#endif
#if 0
  printf("------------------\nmalloc(A) (A cells):\n");
  for (int j=0; j<N; j++) 
    ap = (A *)malloc(sizeof(A));
#endif
  report(N,Asize);
}

//__INSERT_LICENSE__
// $Id: trymem.cpp,v 1.1 2004/12/25 21:12:14 mstorti Exp $

#include <cstdio>
#include <cassert>
#include <sys/resource.h>
#include <sys/vtimes.h>
#include <sys/unistd.h>
#include <cstdlib>
#include <vector>

struct A {
  int k;
  double v;
};

int main() {
  int N = 10000000;
  int Asize = sizeof(A);
  vector<A> va;
  for (int j=0; j<N; j++) new A;
  printf("%d x %dbytes objects generated. OK?\n",
	 N,Asize);
  char file[1000], *line;
  int  mem;
  size_t n=0;
  sprintf(file,"/proc/%d/status",getpid());
  FILE *fid = fopen(file,"r");
  if (0) {
    while(1) {
      assert(getline(&line,&n,fid)!=-1);
      if (sscanf(line,"VmRSS: %d kB",&mem)) break;
    }
  } else {
    va.resize(N);
  }
  fclose(fid);
  double Aused = double(mem)*1024.0/N;
  printf("%dK used, %f bytes/struct approx, excedent %f\n",mem,
	 Aused,Aused-Asize);
}

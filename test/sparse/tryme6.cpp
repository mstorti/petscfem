/*__INSERT_LICENSE__*/
// $Id: tryme6.cpp,v 1.4 2002/07/21 19:01:48 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <cmath>
#include <deque>
#include <vector>
#include <set>

#include "graphdv.h"
#include "graphs.h"

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

void v_print(dvector<int> v) {
  for (int j=0; j<v.size(); j++) 
    printf("%d\n",v.ref(j));
}

int main(int argc, char **argv) {

  graphdv g;
  graph_stl gg;
  vector<int> ngbrs,nngbrs;
  int kk;
  int N = 100000, M = int(N/10), NN = int(sqrt(N/2));
  double fill = double(N)/double(NN*NN)*100.;
  printf("inserting N=%d pairs, in range 1 to NN=%d in chunks of M=%d."
	 " Raw fill %f7.1\%\n",N,NN,M,fill);
  for (kk=1; kk<=N; kk++) { 
    int k = irand(1,NN);
    int l = irand(1,NN);
    // printf("inserting (%d,%d)\n",k,l);
    g.add(k,l); 
    gg.add(k,l); 
    if (kk % M == 0 ) {
      if (g.size()!=gg.size()) {
	printf("on kk=%d bad: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
	g.print("g: ");
	gg.print("gg: ");
	exit(0);
      } else {
	fill = double(gg.size())/double(NN*NN)*100.;
	printf("on kk=%d OK: g.size(): %d, gg.size() %d, fill %5.1f\%\n",
	       kk, g.size(),gg.size(),fill);
      }
    }
  }

  // Check neighbors
  printf("checking neighbors...\n");
  for (int k=1; k<=NN; k++) {
    ngbrs.clear();
    nngbrs.clear();
    g.set_ngbrs(k,ngbrs);
    gg.set_ngbrs(k,nngbrs);
#if 0
    printf("ngbrs of %d: ",k);
    for (vector<int>::iterator q=ngbrs.begin(); q!=ngbrs.end(); q++) 
      printf("%d ",*q);
    printf("\n");
#endif
    
    // Check by comparison with graph_stl
    assert(ngbrs.size()==nngbrs.size());
    for (int q=0; q<ngbrs.size(); q++)
      assert(ngbrs[q] == nngbrs[q]);
  }
  printf("Done.\n");

  g.clear();
  gg.clear();
  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
}

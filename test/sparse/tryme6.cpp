/*__INSERT_LICENSE__*/
// $Id: tryme6.cpp,v 1.3 2002/07/21 16:00:36 mstorti Exp $

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
  int kk;
  int N = 100000, M = int(N/10), NN = int(sqrt(N/2));
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
      } else printf("on kk=%d OK: g.size(): %d, gg.size() %d\n",kk, g.size(),gg.size());
    }
  }
  g.clear();
  gg.clear();
  printf("g.size(): %d, gg.size() %d\n",g.size(),gg.size());
}

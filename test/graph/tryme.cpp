/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.5 2001/11/24 01:01:52 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>

/// Cyclic `rem' 
int crem(int j, int m) {
  int mm = (m>0 ? m : -m);
  if (j>=0 ) {
    return j % mm;
  } else {
    return mm + (j % mm);
  }
}

class TGraph : public Graph {
public:
  int N;
  void set_ngbrs(int elem,vector<int> &ngbrs_v);
  ~TGraph() {clear();}
};

void TGraph::set_ngbrs(int elem,vector<int> &ngbrs_v) {
  ngbrs_v.push_back(crem(elem+1,N));
  ngbrs_v.push_back(crem(elem-1,N));
}

int main(int argc, char **args) {
  const int N=200;
  int k=0;
  while (1) {
    printf("partitioning graph %d\n",k++);
    TGraph *GG = new TGraph;
    Graph &G = *GG;
    GG->N = N;
    G.part(N,N,2);
    delete GG;
//      for (int j=0; j<N; j++) 
//        printf("vrtx %d in proc %d\n",j,G.vrtx_part(j));
  }
}

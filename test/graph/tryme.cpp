/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.3 2001/11/23 20:53:07 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>

class TGraph : public Graph {
public:
  int N;
  void set_ngbrs(int elem,vector<int> &ngbrs_v);
  ~TGraph() {clear();}
};

void TGraph::set_ngbrs(int elem,vector<int> &ngbrs_v) {
  ngbrs_v.push_back(elem+1 % N);
  ngbrs_v.push_back(elem-1 % N);
}

int main(int argc, char **args) {
  const int N=20;
  TGraph GG; 
  Graph &G = GG;
  GG.N = N;
  G.part(N,N,2);
}

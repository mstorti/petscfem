/*__INSERT_LICENSE__*/
// $Id: tryme.cpp,v 1.2 2001/11/23 02:12:23 mstorti Exp $

#include <src/utils.h>
#include <src/graph.h>

class TGraph : public Graph {
  
  void set_ngbrs(int elem,vector<int> &ngbrs_v);
}

void TGraph::set_ngbrs(int elem,vector<int> &ngbrs_v) {
  ngbrs_v.push_back(elem+1 % nelemfat);
  ngbrs_v.push_back(elem-1 % nelemfat);
}

int main(int argc, char **args) {
  const int N=20;
  Graph G(N,N); 
  G.part(2);
}

/*__INSERT_LICENSE__*/
// $Id: tryme3.cpp,v 1.1 2002/07/23 10:33:38 mstorti Exp $

#include <src/utils.h>
#include <src/linkgraph.h>

int main(int argc, char **args) {
  LinkGraphDis graph;
  int M=10;
  graph.init(M);
  for (j=0; j<M; j++) {
    graph.add(j,modulo(j+1,M));
    graph.add(j,modulo(j-1,M));
  }
}

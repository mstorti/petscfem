/*__INSERT_LICENSE__*/
// $Id: tryme5.cpp,v 1.1 2002/07/19 20:16:22 mstorti Exp $

#include <cassert>
#include <cstdio>
#include <cmath>
#include <iostream>                  // for std::cout
#include <utility>                   // for std::pair
#include <algorithm>                 // for std::for_each
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

double drand() {  
  return ((double)(rand()))/((double)(RAND_MAX));
}

int irand(int imin,int imax) {
  return int(rint(drand()*double(imax-imin+1)-0.5))+imin;
}

using namespace boost;

int main(int,char*[]) {
  // create a typedef for the Graph type
  typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

  const int M = 1000;
  // declare a graph object
  Graph g(M);

  // add the edges to the graph object
  for (int i = 0; i < 20*M; ++i)
    add_edge(irand(1,M)-1, irand(1,M)-1, g);

  graph_traits<Graph>::adjacency_iterator ai, ai_end;
  ai = adjacent_vertices(0, g).first();
  ai_end = adjacent_vertices(0, g).first().second;
  while (ai!=ai_end) printf("%d ",*ai++);
  printf("\n");
  return 0;
}
